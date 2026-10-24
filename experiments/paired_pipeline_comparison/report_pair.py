"""Compare completed frozen arms using the repository's existing reporting tool.

This adds provenance, per-stage wall time and validation-status accounting. It
does not add a calling rule or interpret additional assignments as true positives.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import subprocess
import time
import tomllib
from pathlib import Path

p=argparse.ArgumentParser()
for name in ('root','source','retained'):p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--sample',required=True)
a=p.parse_args();d=a.root/'comparisons'/a.sample
complete=json.loads((d/'complete.json').read_text())
assert complete['both_arms_complete'] and complete['downstream_parameters_equal']
gate=json.loads((d/'synthetic-control-gate.json').read_text());assert gate['passed']
out=d/'comparison-report';out.mkdir(exist_ok=False)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def table(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:
        return list(csv.DictReader(f,delimiter='\t'))

project=tomllib.loads((a.source/'pyproject.toml').read_text())
expected_files=project.get('tool',{}).get('mgs-workflow',project)['expected-outputs-downstream']
assert expected_files
missing=[]
for arm in ('production','protein_variant'):
    for relative in expected_files:
        path=d/(arm+'-downstream/output')/relative.replace('{GROUP}',a.sample)
        if not path.is_file():missing.append(str(path))
save(out/'expected-files.json',{'missing':missing,'expected_per_arm':len(expected_files)})
if missing:raise RuntimeError('Missing declared downstream output: '+str(missing))

# The existing reporter expects INDEX/output/results, whereas the retained
# reference is mounted at reference/results. This view changes only the path.
view=out/'reference-view/output';view.mkdir(parents=True)
(view/'results').symlink_to(a.retained/'reference/results',target_is_directory=True)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/pandas:91e85c3fd8e9026c'
image_id=subprocess.check_output(['docker','image','inspect',image,'--format','{{.Id}}'],text=True).strip()
cmd=['docker','run','--rm','--network','none','-v',f'{a.root}:{a.root}',
     '-v',f'{a.retained}:{a.retained}:ro','--entrypoint','python',image_id,
     str(a.source/'bin/compare_downstream_runs.py'),'--reference',str(d/'production-downstream/output'),
     '--candidate',str(d/'protein_variant-downstream/output'),
     '--candidate-index',str(view.parent),'--reference-index',str(view.parent),
     '--out',str(out/'standard'),'--schema-dir',str(a.source/'schemas'),
     '--pyproject',str(a.source/'pyproject.toml')]
started=datetime.datetime.now(datetime.timezone.utc).isoformat();clock=time.monotonic()
with (out/'standard-report.log').open('w') as f:result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
save(out/'execution.json',{'command':cmd,'started_utc':started,'seconds':time.monotonic()-clock,
    'exit_code':result.returncode,'image_id':image_id,'pipeline_source_commit':complete['source_commit'],
    'report_source_hashes':{name:sha(a.source/'bin'/name) for name in ('compare_downstream_runs.py','downstream_metrics.py')}})
if result.returncode:raise RuntimeError('Existing comparison report failed; retained log and receipt')
missing_columns=[r for r in table(out/'standard/column_conformance.tsv')
                 if r['missing_vs_schema_reference'] or r['missing_vs_schema_candidate']]
skipped=table(out/'standard/skipped_groups.tsv')
save(out/'integrity.json',{'passed':not missing_columns and not skipped,'missing_columns':missing_columns,'skipped_groups':skipped})
if missing_columns or skipped:raise RuntimeError('Incomplete comparison inventory; see retained integrity.json')

stages=[]
for path in sorted(d.glob('*.json')):
    r=json.loads(path.read_text())
    if 'stage' in r and 'finished_utc' in r:
        assert r['exit_code']==0,path
        stages.append({k:r[k] for k in ('stage','started_utc','finished_utc','wall_seconds')})
stages.sort(key=lambda r:r['started_utc'])
with (out/'stage-wall-times.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['stage','started_utc','finished_utc','wall_seconds'],delimiter='\t');w.writeheader();w.writerows(stages)
statuses={};identities={};inventory=[]
for arm in ('production','protein_variant'):
    result_dir=d/(arm+'-downstream/output/results_downstream')
    files=list(result_dir.glob('*_validation_hits.tsv.gz'));assert len(files)==1,files
    rows=table(files[0]);ids=[(r['group'],r['seq_id']) for r in rows]
    assert len(set(ids))==len(ids)
    identities[arm]=set(ids)
    statuses[arm]={'rows':len(rows),'validation_status':dict(collections.Counter(r['validation_status'] for r in rows)),
        'joint_exemplars':sum(r['seq_id']==r['prim_align_dup_exemplar']==r['sim_dup_exemplar'] for r in rows)}
    for f in sorted(result_dir.iterdir()):
        if f.is_file():inventory.append({'arm':arm,'file':str(f),'bytes':f.stat().st_size,'sha256':sha(f)})
save(out/'validation-status-counts.json',statuses)
save(out/'output-checksums.json',inventory)
save(out/'accounting.json',{'completed':True,'sample':a.sample,'pairs_per_arm':complete['pairs'],
    'shared_output_pair_ids':len(identities['production']&identities['protein_variant']),
    'production_only_output_pair_ids':len(identities['production']-identities['protein_variant']),
    'protein_only_output_pair_ids':len(identities['protein_variant']-identities['production']),
    'sum_stage_wall_seconds':sum(r['wall_seconds'] for r in stages),
    'cost_interpretation':'Stage wall time excludes full-source acquisition, boot, checkpoints and other worker time. Billable worker time must be accounted separately. No cost forecast is inferred from a single case.',
    'classification_interpretation':'The standard pipeline outputs are compared without an added NR agreement or three-exemplar calling rule. Changes are not ground-truth sensitivity or specificity.'})
print(json.dumps({'report':str(out),'validation_status_counts':statuses},indent=2))

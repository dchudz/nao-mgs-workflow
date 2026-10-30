"""Run original validation on pooled, independently duplicate-marked sample tables."""
import csv,datetime,gzip,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');s=r/'cohort-storage';d=s/'validation-cohort';d.mkdir(exist_ok=True)
planned=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())]
cutoff=datetime.datetime.fromisoformat('2026-10-23T05:00:00+00:00')
while True:
 found={a:p.parent for p in (s/'batches').glob('*/complete.json') for a in json.loads(p.read_text())['samples']}
 if all(a in found for a in planned):break
 if datetime.datetime.now(datetime.timezone.utc)>=cutoff:raise RuntimeError('Full 1000-sample validation cannot start within the budget safeguard; no partial cohort substituted.')
 time.sleep(30)
samples=[a for a in planned if a in found];assert samples
scope={'planned':len(planned),'completed_samples':samples,'missing_samples':[a for a in planned if a not in found],'duplicate_unit':'Biological sample and arm, before pooling.','validation_unit':'20 selected exemplars per species per pooled arm, using the existing workflow.','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(d/'scope.json').write_text(json.dumps(scope,indent=2))
with (d/'manifest.tsv').open('w') as manifest:
 manifest.write('id\tpath\n')
 for arm in ['nt','nt22','nt20','aa12','union']:
  path=d/f'{arm}-pooled.tsv.gz';header=None
  with gzip.open(path,'wb') as out:
   for a in samples:
    files=list((found[a]/'output/marked').glob(f'*_{a}_{arm}_duplicate_reads.tsv.gz'));assert len(files)==1,(a,arm)
    with gzip.open(files[0],'rb') as source:
     first=source.readline()
     if header is None:header=first;out.write(first)
     else:assert header==first
     for chunk in iter(lambda:source.read(1024*1024),b''):out.write(chunk)
  manifest.write(f'{arm}\t{path}\n')
cmd=[str(r/'nextflow'),'-log',str(d/'nextflow.log'),'run','experiments/amino_acids/validation-pooled.nf','-c','experiments/amino_acids/validation-cohort.config','-work-dir',str(d/'work'),'--manifest',str(d/'manifest.tsv'),'--ref_dir',str(r/'reference'),'--output',str(d/'output')]
remaining=(datetime.datetime.fromisoformat('2026-10-23T14:30:00+00:00')-datetime.datetime.now(datetime.timezone.utc)).total_seconds()
assert remaining>600,'Insufficient validation time before cleanup deadline'
with (d/'run.log').open('w') as log:subprocess.run(cmd,cwd=r,stdout=log,stderr=log,check=True,timeout=min(28800,remaining))
assert len(list((d/'output/COPY_HITS').glob('*_validation_hits.tsv.gz')))==5
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-cohort/','--recursive','--exclude','work/*','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_samples':len(samples),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),'s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-cohort/complete.json','--only-show-errors'],check=True)
print('COHORT_VALIDATION_COMPLETE',flush=True)

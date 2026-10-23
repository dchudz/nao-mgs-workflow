"""Audit closed output from a time-limited process without changing its status."""
import argparse, collections, csv, datetime, hashlib, json, re, shutil, subprocess
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--commit',required=True);p.add_argument('--output-name',default='output-audit-v2');a=p.parse_args()
assert re.fullmatch('[0-9a-f]{40}',a.commit)
assert re.fullmatch('output-audit(?:-v[0-9]+)?',a.output_name)
r=Path('/home/ec2-user/experiment');protein=r/'protein'
d=protein/'seed-nr-forced-double-index-trace';out=d/a.output_name
resource=json.loads((d/'resources.json').read_text())
protocol=json.loads((d/'search-protocol.json').read_text())
assert resource['status']=='time_limit' and resource['returncode']==-15
assert not (d/'complete.json').exists()
assert not subprocess.run(['pgrep','-x','diamond'],capture_output=True).stdout.strip()
log=(d/'nr-search.log').read_text()
assert 'Algorithm: Double-indexed' in log and 'Closing the output file...  [0.001s]' in log
assert {int(x) for x in re.findall(r'reference block (\d+)/221',log)}==set(range(1,222))
reported=re.search(r'Reported (\d+) pairwise alignments, (\d+) HSPs\.',log)
assert reported and reported[1]==reported[2]
queries=[x[1:].split()[0] for x in (d/'queries.fna').read_text().splitlines() if x.startswith('>')]
assert len(queries)==len(set(queries))==1436
query_set=set(queries)
assert hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest()==protocol['query_sha256']=='b5d74d481f0291a671af931ff0e14f13b97f58f38ecd85b2bd2950e6cb3b841a'
out.mkdir(exist_ok=False)
for name in ('queries.fna','seed-selection.tsv'):shutil.copyfile(d/name,out/name)
raw=d/'nr-with-ordinals.tsv';initial=raw.stat();h=hashlib.sha256();projected=hashlib.sha256();rows=0;new=collections.Counter();human=[]
wanted={'SEED_SRR33038061.678928042/'+str(m) for m in (1,2)}|{'SEED_SRR35987653.689121697/'+str(m) for m in (1,2)}
with raw.open('rb') as inp,(out/'nr.tsv').open('wb') as dest:
    for line in inp:
        h.update(line);f=line.rstrip(b'\n').split(b'\t');assert len(f)==8 and f[7].isdigit()
        q=f[0].decode();assert q in query_set
        row=b'\t'.join(f[:7])+b'\n';dest.write(row);projected.update(row);new[row.rstrip(b'\n').decode()]+=1;rows+=1
        if q in wanted and b'9606' in f[6].split(b';'):
            human.append(dict(zip('query subject bitscore evalue alignment_AA identity_pct taxids global_ordinal'.split(),[x.decode() for x in f])))
assert raw.stat().st_size==initial.st_size and raw.stat().st_mtime_ns==initial.st_mtime_ns
assert rows==int(reported[1])
projection=dict(rows=rows,raw_bytes=initial.st_size,raw_sha256=h.hexdigest(),projected_sha256=projected.hexdigest(),raw_fields=8,analysis_fields=7,logged_output_closed=True,logged_reference_blocks=221)
(out/'projection.json').write_text(json.dumps(projection,indent=2)+'\n')
comparisons={}
for label,path in [('pooled_double',protein/'combined-validation/seed-alignments.tsv'),('standalone_auto_query',protein/'seed-nr-standalone/nr.tsv')]:
    old=collections.Counter();original_hash=hashlib.sha256()
    with path.open('rb') as f:
        for line in f:
            original_hash.update(line);old[line.rstrip(b'\r\n').decode()]+=1
    left=old-new;right=new-old
    comparisons[label]=dict(alignment_multiset_equal=old==new,original_sha256=original_hash.hexdigest(),original_rows=sum(old.values()),audited_rows=rows,rows_only_original=sum(left.values()),rows_only_audited=sum(right.values()),raw_differing_query_mates=len({x.split('\t',1)[0] for x in left.keys()|right.keys()}))
    del old,left,right
(out/'context-comparison.json').write_text(json.dumps(comparisons,indent=2)+'\n')
assert comparisons['pooled_double']['alignment_multiset_equal'], 'Retain raw audit, but do not classify potentially incomplete output as a completed search'
# Exact raw-output equality permits replaying the existing classifier as an artifact
# comparison. It does not turn the interrupted process into a successful search.
subprocess.run(['python3',str(r/'summarize-seed-nr.py'),str(out)],check=True)
# The existing stage summarizer requires clean GNU exit statistics. Preserve
# that guard; this interrupted process has no such timing record.
def calls(path):
    values=list(csv.DictReader(path.open(),delimiter='\t'));result={x['query_pair']:x for x in values};assert len(values)==len(result)==746;return result
after=calls(out/'seed-validation.tsv');fields=['nr_status','nr_taxid_lca','aligned_mates','taxonomic_relation'];changes=[]
for label,path in [('pooled_double',protein/'combined-validation/seed-validation.tsv'),('standalone_auto_query',protein/'seed-nr-standalone/seed-validation.tsv')]:
    before=calls(path);assert before.keys()==after.keys()
    changed=[q for q in before if any(before[q][k]!=after[q][k] for k in fields)]
    comparisons[label].update(pairs_with_summary_difference=len(changed),status_disagreements=sum(before[q]['nr_status']!=after[q]['nr_status'] for q in before),taxid_disagreements=sum(before[q]['nr_taxid_lca']!=after[q]['nr_taxid_lca'] for q in before))
    changes.extend(dict(comparison=label,query_pair=q,**{'original_'+k:before[q][k] for k in fields},**{'audited_'+k:after[q][k] for k in fields}) for q in changed)
for name,data,fieldnames in [('human-target-ordinals.tsv',human,'query subject bitscore evalue alignment_AA identity_pct taxids global_ordinal'.split()),('context-disagreements.tsv',changes,['comparison','query_pair']+['original_'+k for k in fields]+['audited_'+k for k in fields])]:
    with (out/name).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames,delimiter='\t');writer.writeheader();writer.writerows(data)
(out/'context-comparison.json').write_text(json.dumps(comparisons,indent=2)+'\n')
receipt=dict(audited_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audit_source_commit=a.commit,audit_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),search_source_commit=protocol['source_commit'],search_status_unchanged=resource,projection=projection,comparisons=comparisons,human_hsp_rows=len(human),queries=1436,pairs=746,allocated_time_limit_USD=resource['wall_seconds']/3600*2.679511111111111,logged_search_total_seconds=float(re.search(r'Total time = ([\d.]+)s',log)[1]),limits=['This is a completed artifact audit, not a successful process-completion receipt. The original search remains time_limit, returncode -15.','The log records all 221 blocks, closed output, a final reported row count, and total search time. No clean process exit or final GNU timing statistics were recorded.','All seven original fields and duplicate multiplicity exactly reproduce the earlier pooled diagnostic output. Classifications are replayed only after that equality check.','The fixed-mode standalone artifact removes the observed discrepancies on these 746 selected pairs. This does not establish parity on other inputs or biological accuracy.'])
(out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)

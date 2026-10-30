from pathlib import Path
import json,csv,gzip,collections,hashlib,itertools,subprocess
r=Path('/home/ec2-user/experiment');d=r/'e2e/SRR37076045'
for sub in ['protein-run/published','protein-run/work/62/f6a71d0c492f2c3e0e47e2f547717e']:
 p=d/sub
 for q in p.rglob('*'):
  if q.is_file() and not q.is_symlink() and (q.suffix=='.json' or q.name in ['.command.err','.command.log']):
   print(str(q.relative_to(d)),q.read_text(errors='replace')[-3000:])
left=d/'production-downstream/output/results_downstream/SRR37076045_validation_hits.tsv.gz'
right=d/'validation_control-downstream/output/results_downstream/SRR37076045_validation_hits.tsv.gz'
metrics={}; counters=collections.Counter();rows=0
with gzip.open(left,'rt') as lf,gzip.open(right,'rt') as rf:
 l=csv.DictReader(lf,delimiter='\t');rr=csv.DictReader(rf,delimiter='\t');keys=[k for k in l.fieldnames if not k.startswith('validation_')]
 for x,y in itertools.zip_longest(l,rr):
  assert x is not None and y is not None and all(x[k]==y[k] for k in keys)
  rows+=1
  for arm,z in [('production',x),('validation_control',y)]:
   metrics.setdefault(arm,collections.Counter())[z['validation_status']]+=1
  if x['validation_status']!=y['validation_status']:counters['status_changed']+=1
  if x['validation_staxid_lca']!=y['validation_staxid_lca']:counters['natural_lca_changed']+=1
record={'sample':d.name,'nonvalidation_columns_identical':True,'rows':rows,'statuses':metrics,'changes':counters,'note':'Alignment outcomes, not confirmed virus or infection counts.'}
(d/'native-validation-comparison.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
for arm in ['production','validation_control']:
 trace=d/(arm+'-downstream')/'trace.tsv'
 for row in csv.DictReader(trace.open(),delimiter='\t'):
  if row.get('name','').endswith('BLAST (1)') or ':BLAST ' in row.get('name',''):print(arm,'BLAST',row)

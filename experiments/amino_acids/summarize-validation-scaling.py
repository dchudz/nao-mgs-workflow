"""Report search work and whole-worker allocation without claiming incremental billing."""
import argparse,collections,gzip,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root');a=p.parse_args();d=Path(a.root)
rows=json.loads((d/'results.json').read_text())
# Linux/shared x2iedn.4xlarge plus its configured gp3 allocation, using recorded AWS rates.
compute_hour=3.3345;storage_hour=(2048*.08+3000*.005+875*.04)/730
for x in rows:
 text=(d/f"run{x['run']}-cap{x['cap']}"/'blast.time').read_text()
 def number(label):
  m=re.search(re.escape(label)+r':\s*([\d.]+)',text);assert m,label;return float(m[1])
 x['cpu_seconds']=number('User time (seconds)')+number('System time (seconds)')
 x['peak_rss_kib']=int(number('Maximum resident set size (kbytes)'))
 x['allocated_worker_cost_usd']=x['wall_seconds']/3600*(compute_hour+storage_hour)
def fasta(path):
 result={};q=None
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   q=line[1:].split()[0];assert q not in result;result[q]=''
  else:
   assert q is not None;result[q]+=line
 return result
queries={cap:fasta(d/f'cap{cap}/queries.fna') for cap in [1,5,20,100]}
parity=[]
for small,large in [(1,5),(5,20),(20,100)]:
 assert queries[small].keys()<=queries[large].keys()
 changed=[q for q,seq in queries[small].items() if queries[large][q]!=seq]
 parity.append({'smaller_cap':small,'larger_cap':large,'shared_queries':len(queries[small]),'changed_sequences':len(changed),'changed_query_ids':changed})
replicas=[x for x in rows if x['cap']==20]
assert len(replicas)==2
observed=[]
for x in replicas:
 with gzip.open(d/f"run{x['run']}-cap20"/'hits.tsv.gz','rt') as f:observed.append(collections.Counter(f))
repeat={'cap':20,'first_rows':sum(observed[0].values()),'second_rows':sum(observed[1].values()),'alignment_multiset_equal':observed[0]==observed[1],'rows_only_first':sum((observed[0]-observed[1]).values()),'rows_only_second':sum((observed[1]-observed[0]).values())}
summary={'runs':rows,'full_worker_hour_usd':compute_hour+storage_hour,'caps':[1,5,20,100],'shared_query_sequence_checks':parity,'repeated_cap20_alignment_check':repeat,'notes':['Allocated cost charges the full instance and EBS for measured wall time; concurrent jobs share that allocation, so do not add it to their attributed cost.','Warm shared reference cache, concurrent native tasks, and non-isolated wall time. The repeated cap20 run measures temporal variability.','No query-count interpolation, no extrapolated verdicts, and no claim that unaligned reads are negatives.']}
(d/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='runs'}));print(json.dumps([{k:v for k,v in x.items() if k!='command'} for x in rows]))

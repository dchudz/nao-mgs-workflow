"""Measure original-read membership in frozen nucleotide screens for every protein candidate."""
import collections,csv,datetime,gzip,hashlib,json,subprocess,sys,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
source=r/'protein-repeat-gate';d=r/'protein-screen-boundaries';d.mkdir(exist_ok=True)
assert not (d/'complete.json').exists(),'Do not overwrite a completed audit'
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
stop=datetime.datetime.fromisoformat('2026-10-23T13:30:00+00:00')
while not (source/'complete.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Budget cleanup reserve reached before input ready'
 time.sleep(30)
viral={x['query_pair'] for x in csv.DictReader((source/'viral-pairs.tsv').open(),delimiter='\t')};members=collections.defaultdict(list);member_to_query={};by_sample=collections.defaultdict(set)
with gzip.open(source/'all.members.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  q=x['query_pair'];sid=x['seq_id']
  if q not in viral:continue
  assert sid not in member_to_query;member_to_query[sid]=q;members[q].append(sid);by_sample[q.split('_',1)[0]].add(sid)
assert members.keys()==viral
order=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200];assert set(by_sample)<=set(order)
retained={arm:set() for arm in ['nt','nt20']};scan_counts=[];start=time.monotonic()
for sample in order:
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Budget cleanup reserve reached during audit'
 for arm in ['nt','nt20']:
  raw=d/'incoming.fastq';subprocess.run(['aws','--profile','pipeline-codex','s3','cp',uri+f'cohort-screens/{sample}/{arm}.fastq',str(raw),'--only-show-errors'],check=True)
  n=0
  for key,record in pairs(raw):
   n+=1;sid=key[0].decode().lstrip('@')
   if sid in by_sample[sample]:assert sid not in retained[arm];retained[arm].add(sid)
  scan_counts.append(dict(sample=sample,arm=arm,screened_candidate_pairs_scanned=n));raw.unlink()
 print(json.dumps({'audited_samples':len(scan_counts)//2,'planned':200}),flush=True)
rows=[]
for q,ids in sorted(members.items()):
 x=dict(query_pair=q,sample=q.split('_',1)[0],original_members=len(ids))
 for arm in ['nt','nt20']:
  n=sum(s in retained[arm] for s in ids);x[arm+'_retained_original_members']=n;x[arm+'_screen_status']='none_retained' if n==0 else 'all_retained' if n==len(ids) else 'some_retained'
 rows.append(x)
with (d/'query-membership.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with (d/'original-member-screen-membership.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['query_pair','seq_id','present_in_NT24_screen','present_in_NT20_screen'],delimiter='\t');w.writeheader()
 for sid,q in sorted(member_to_query.items()):w.writerow(dict(query_pair=q,seq_id=sid,present_in_NT24_screen=sid in retained['nt'],present_in_NT20_screen=sid in retained['nt20']))
summary={'samples':200,'candidate_queries':len(rows),'original_members':len(member_to_query),'wall_seconds':time.monotonic()-start,'arms':{arm:{'retained_original_members':len(retained[arm]),'query_screen_status':dict(collections.Counter(x[arm+'_screen_status'] for x in rows))} for arm in retained},'interpretation':'Read-ID membership in saved raw nucleotide-screen FASTQs, not a new alignment or validation verdict. Final native assignment and NR status are joined separately. Exact-sequence query membership is not independent validation of every duplicate original member.'}
(d/'summary.json').write_text(json.dumps(summary,indent=2));(d/'scan-counts.json').write_text(json.dumps(scan_counts,indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri+'protein-screen-boundaries/','--recursive','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),uri+'protein-screen-boundaries/complete.json','--only-show-errors'],check=True);print(json.dumps(summary),flush=True)

"""Post-hoc raw-gate attrition check; never alter the frozen protein calls."""
import collections,csv,datetime,gzip,hashlib,json,subprocess,sys,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'))
from union_pairs import pairs
d=r/'protein-repeat-gate';d.mkdir(exist_ok=True)
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
deadline=datetime.datetime.fromisoformat('2026-10-23T05:00:00+00:00')
def cp(source,destination):subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(source),str(destination),'--only-show-errors'],check=True)
while True:
 check=subprocess.run(['aws','--profile','pipeline-codex','s3api','head-object','--bucket','pipeline-lab-839993661825-us-east-1','--key','runs/amino-acid-20261021/producer-finished.json'],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
 if check.returncode==0:break
 assert '404' in check.stderr or 'Not Found' in check.stderr,check.stderr
 if datetime.datetime.now(datetime.timezone.utc)>deadline:raise RuntimeError('Primary-screen completion deadline reached; no sample substitution')
 time.sleep(30)
cp(uri+'protein200/viral-pairs.tsv',d/'viral-pairs.tsv')
cp(uri+'protein200/all.members.tsv.gz',d/'all.members.tsv.gz')
cp(uri+'assets/virus.aa12-repeat-filtered.bin',d/'repeat-index.bin')
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()
assert digest(d/'repeat-index.bin')=='e4c86f4547c92f6d73dd45bbbd7e1b923c6a2dd28bda7553905ca011bd0febb9'
calls={x['query_pair']:x for x in csv.DictReader(open(d/'viral-pairs.tsv'),delimiter='\t')}
member_to_query={};members=collections.defaultdict(set);by_sample=collections.defaultdict(set)
with gzip.open(d/'all.members.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  q=x['query_pair'];sid=x['seq_id']
  if q not in calls:continue
  assert sid not in member_to_query;member_to_query[sid]=q;members[q].add(sid);by_sample[q.split('_',1)[0]].add(sid)
assert set(members)==set(calls)
order=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200]
assert set(by_sample)<=set(order)
counts=[]
with (d/'candidate-members.fastq').open('xb') as selected:
 for sample in order:
  wanted=by_sample[sample]
  if not wanted:counts.append({'sample':sample,'candidate_members':0});continue
  raw=d/'incoming.fastq';cp(uri+f'cohort-screens/{sample}/aa12.fastq',raw);seen=set()
  for key,record in pairs(raw):
   sid=key[0].decode().lstrip('@')
   if sid in wanted:
    assert sid not in seen;seen.add(sid);selected.write(record)
  assert seen==wanted,(sample,len(wanted-seen));raw.unlink();counts.append({'sample':sample,'candidate_members':len(seen)})
  print(json.dumps({'extracted_samples':len(counts),'planned':200}),flush=True)
cmd=[str(r/'rust-tools/target-speed/release/screen_amino_acids'),'screen','--index',str(d/'repeat-index.bin'),'--reads',str(d/'candidate-members.fastq'),'--matched',str(d/'retained-members.fastq'),'--threads','8']
with (d/'screen.log').open('w') as log:subprocess.run(['/usr/bin/time','-v']+cmd,stdout=log,stderr=log,check=True,timeout=3600)
retained={key[0].decode().lstrip('@') for key,raw in pairs(d/'retained-members.fastq')};assert retained<=member_to_query.keys()
out=[]
for q,ids in members.items():
 n=len(ids&retained);out.append(dict(query_pair=q,sample=q.split('_',1)[0],original_members=len(ids),raw_members_retained=n,raw_members_excluded=len(ids)-n,gate_status='none_retained' if n==0 else 'all_retained' if n==len(ids) else 'some_retained'))
with (d/'query-gate.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary={'samples':200,'candidate_exact_pairs':len(out),'original_candidate_members':len(member_to_query),'retained_original_members':len(retained),'gate_status':dict(collections.Counter(x['gate_status'] for x in out)),'candidate_members_sha256':digest(d/'candidate-members.fastq'),'retained_members_sha256':digest(d/'retained-members.fastq'),'command':cmd,'interpretation':'Post-hoc raw-screen attrition diagnostic on original members of the fixed first-200 protein candidates. This does not rerun FASTP, contaminant filters or alignments and does not establish final-pipeline parity. A query retains raw support if at least one original member passes the repeat-filtered index. Original NR outcomes are joined later without changing selection.'}
(d/'summary.json').write_text(json.dumps(summary,indent=2));(d/'sample-counts.json').write_text(json.dumps(counts,indent=2));(d/'repeat-index.bin').unlink()
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri+'protein-repeat-gate/','--recursive','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));cp(d/'complete.json',uri+'protein-repeat-gate/complete.json');print(json.dumps(summary),flush=True)

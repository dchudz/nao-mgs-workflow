"""Measure the frozen diagnostic alone after primary search, backup and analysis complete."""
import collections,csv,datetime,hashlib,json,os,signal,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');primary=r/'protein/combined-validation';d=r/'protein/seed-nr-standalone';d.mkdir(exist_ok=True)
assert not (d/'nr.tsv').exists(),'Do not overwrite a started search'
source=primary/'combined-queries.fna';assert hashlib.sha256(source.read_bytes()).hexdigest()=='6aa28bedbaf612c3a7c5a895b24a3be4f18d2c8264a439846427d943b289c597'
selected=[];record=[]
with source.open('rb') as f:
 for line in f:
  if line.startswith(b'>'):
   if record and record[0].startswith(b'>SEED_'):selected.extend(record)
   record=[line]
  else:record.append(line)
 if record and record[0].startswith(b'>SEED_'):selected.extend(record)
query=d/'queries.fna';query.write_bytes(b''.join(selected));names=[line[1:].split()[0].decode() for line in selected if line.startswith(b'>')];assert len(names)==len(set(names))==1436
selection=primary/'seed-selection.tsv';rows=list(csv.DictReader(selection.open(),delimiter='\t'));assert len(rows)==746 and sum(int(x['available_mates']) for x in rows)==1436
assert {q.rsplit('/',1)[0].removeprefix('SEED_') for q in names}=={x['seq_id'] for x in rows}
(d/'seed-selection.tsv').write_bytes(selection.read_bytes());addendum=json.loads((r/'standalone-seed-nr-addendum.json').read_text());(d/'protocol.json').write_text(json.dumps(dict(addendum,queries=len(names),query_sha256=hashlib.sha256(query.read_bytes()).hexdigest(),prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/seed-nr-standalone/'
def upload():subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri,'--recursive','--exclude','complete.json','--only-show-errors'],check=True)
upload();print('STANDALONE_SEED_NR_PREPARED_WAITING_FOR_PRIMARY_ANALYSIS',flush=True)
budget_stop=datetime.datetime.fromisoformat('2026-10-23T08:30:00+00:00')
while not (r/'nr-primary-analysis-complete.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<budget_stop-datetime.timedelta(minutes=10),'No time remains before budget cleanup reserve'
 time.sleep(30)
assert (primary/'complete.json').exists()
subprocess.run(['aws','--profile','pipeline-codex','s3api','head-object','--bucket','pipeline-lab-839993661825-us-east-1','--key','runs/amino-acid-20261021/protein-nr-combined/complete.json'],check=True,stdout=subprocess.DEVNULL)
cmd=json.loads((primary/'search-protocol.json').read_text())['command'].copy()
assert cmd[cmd.index('--threads')+1]=='32' and cmd[cmd.index('--index-chunks')+1]=='4'
cmd[cmd.index('--query')+1]=str(query);cmd[cmd.index('--out')+1]=str(d/'nr.tsv')
protocol={'command':cmd,'queries':len(names),'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'instance_type':'r7a.8xlarge','memory_GiB':256,'worker_hourly_cost_usd':2.4344+(1280*.08+3000*.005+875*.04)/730,'max_runtime_seconds':(budget_stop-datetime.datetime.now(datetime.timezone.utc)).total_seconds()};assert protocol['max_runtime_seconds']>600
(d/'search-protocol.json').write_text(json.dumps(protocol,indent=2));print(json.dumps(protocol),flush=True)
start=time.monotonic()
with (d/'nr-search.log').open('w') as log:
 proc=subprocess.Popen(['/usr/bin/time','-v']+cmd,stdout=log,stderr=log,start_new_session=True)
 try:
  code=proc.wait(timeout=protocol['max_runtime_seconds'])
  if code:raise subprocess.CalledProcessError(code,cmd)
 except BaseException as error:
  if proc.poll() is None:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=15)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
  (d/'incomplete.json').write_text(json.dumps({'error':repr(error),'wall_seconds':time.monotonic()-start,'primary_results_unchanged':True},indent=2));upload();raise
(d/'resources.json').write_text(json.dumps({'wall_seconds':time.monotonic()-start},indent=2))
subprocess.run(['python3',str(r/'summarize-seed-nr.py'),str(d)],check=True)
subprocess.run(['python3',str(r/'summarize-diamond-stages.py'),str(d/'nr-search.log'),str(d/'stage-timings.json')],check=True)
def alignments(path):
 with path.open() as f:return collections.Counter(line.rstrip('\n') for line in f)
old=alignments(primary/'seed-alignments.tsv');new=alignments(d/'nr.tsv');only_old=old-new;only_new=new-old
raw_queries=sorted({line.split('\t',1)[0] for line in only_old.keys()|only_new.keys()})
def calls(path):
 rows=list(csv.DictReader(path.open(),delimiter='\t'));out={x['query_pair']:x for x in rows};assert len(out)==len(rows);return out
before=calls(primary/'seed-validation.tsv');after=calls(d/'seed-validation.tsv');assert before.keys()==after.keys() and len(after)==746
changes=[]
for q in sorted(before):
 fields=[k for k in ['nr_status','nr_taxid_lca','aligned_mates','taxonomic_relation'] if before[q][k]!=after[q][k]]
 if fields:changes.append(dict(query_pair=q,changed_fields=';'.join(fields),**{'pooled_'+k:before[q][k] for k in fields},**{'standalone_'+k:after[q][k] for k in fields}))
fields=['query_pair','changed_fields']+['pooled_'+k for k in ['nr_status','nr_taxid_lca','aligned_mates','taxonomic_relation']]+['standalone_'+k for k in ['nr_status','nr_taxid_lca','aligned_mates','taxonomic_relation']]
with (d/'context-disagreements.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(changes)
comparison={'pairs':len(after),'queries':len(names),'alignment_multiset_equal':old==new,'pooled_alignment_rows':sum(old.values()),'standalone_alignment_rows':sum(new.values()),'alignment_rows_only_pooled':sum(only_old.values()),'alignment_rows_only_standalone':sum(only_new.values()),'raw_alignment_differing_query_mates':len(raw_queries),'raw_alignment_differing_pair_queries':len({q.rsplit('/',1)[0] for q in raw_queries}),'raw_alignment_differing_query_IDs':raw_queries,'pairs_with_summary_difference':len(changes),'status_disagreements':sum(before[q]['nr_status']!=after[q]['nr_status'] for q in before),'taxid_disagreements':sum(before[q]['nr_taxid_lca']!=after[q]['nr_taxid_lca'] for q in before),'interpretation':'Same query bytes, executable, comprehensive reference and frozen search parameters; only the surrounding query pool changes. Standalone results are secondary and never replace pooled primary results.'}
(d/'context-comparison.json').write_text(json.dumps(comparison,indent=2));upload();(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queries':1436,'pairs':746},indent=2));subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True);print('STANDALONE_SEED_NR_COMPLETE',flush=True)

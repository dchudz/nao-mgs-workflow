"""Foreground, bounded batch processing of atomically published screen results."""
import subprocess,json,time,csv,gzip,datetime,tarfile,shutil,re,os
from pathlib import Path
r=Path('/home/ec2-user/experiment');s=r/'cohort-storage';s.mkdir(exist_ok=True)
rows=json.loads((r/'heldout-ranges.json').read_text());order={x['run_accession']:i for i,x in enumerate(rows)}
bucket='pipeline-lab-839993661825-us-east-1';prefix='runs/amino-acid-20261021/cohort-screens/'
deadline=datetime.datetime.fromisoformat('2026-10-23T04:30:00+00:00')
profile=['--profile','pipeline-codex','--region','us-east-1']
def aws(*args):return subprocess.check_output(['aws',*profile,*args],text=True)
def ready():
 keys=json.loads(aws('s3api','list-objects-v2','--bucket',bucket,'--prefix',prefix,'--query','Contents[?ends_with(Key, `/complete.json`)].Key','--output','json')) or []
 return {k.split('/')[-2] for k in keys}&order.keys()
def completed():
 result=set(json.loads((s/'prior-completed.json').read_text())) if (s/'prior-completed.json').exists() else set()
 for p in (s/'batches').glob('*/complete.json'):result.update(json.loads(p.read_text())['samples'])
 return result
def process_batch(ids,batch):
 d=s/'batches'/f'{os.environ.get("AA_BATCH_PREFIX","s")}{batch:03d}';d.mkdir(parents=True,exist_ok=True);incoming=d/'incoming';incoming.mkdir(exist_ok=True)
 manifest=d/'manifest.tsv';out=d/'output';work=d/'work'
 for a in ids:
  if (incoming/a/'complete.json').exists():continue
  subprocess.run(['aws',*profile,'s3','cp',f's3://{bucket}/{prefix}{a}/',str(incoming/a)+'/', '--recursive','--only-show-errors'],check=True)
 with manifest.open('w') as f:
  f.write('id\tpath\n')
  for a in ids:
   observed=json.loads((incoming/a/'complete.json').read_text());assert observed['observed_pairs']==5000000
   for arm in ['nt','nt22','nt20','aa12','union']:f.write(f'{a}_{arm}\t{incoming/a/(arm+".fastq")}\n')
 start=time.monotonic()
 cmd=[str(r/'nextflow'),'-log',str(d/'nextflow.log'),'run','experiments/amino_acids/heldout-workflow.nf','-c','experiments/amino_acids/heldout.config','-work-dir',str(work),'--manifest',str(manifest),'--ref_dir',str(r/'reference'),'--adapters',str(r/'ref/adapters.fasta'),'--output',str(out)]
 if (d/'nextflow.log').exists():
  ids_found=re.findall(r'Session UUID: ([0-9a-f-]{36})',(d/'nextflow.log').read_text())
  assert len(set(ids_found))==1,ids_found
  cmd.extend(['-resume',ids_found[0]])
  attempt=len(list(d.glob('run-attempt-*.log')))+1
  shutil.copyfile(d/'nextflow.log',d/f'nextflow-attempt-{attempt}.log')
  if (d/'run.log').exists():shutil.copyfile(d/'run.log',d/f'run-attempt-{attempt}.log')
 with (d/'run.log').open('w') as log:subprocess.run(cmd,cwd=r,stdout=log,stderr=log,check=True,timeout=5400)
 # Every completed sample/arm must have the primary result and duplicate-marked result.
 assert len(list((out/'hits').glob('*_virus_hits.tsv.gz')))==len(ids)*5
 assert len(list((out/'marked').glob('*.tsv.gz')))==len(ids)*5
 metrics=[]
 for a in ids:
  for arm in ['nt','nt22','nt20','aa12','union']:
   files=list((out/'marked').glob(f'*_{a}_{arm}_duplicate_reads.tsv.gz'));assert len(files)==1,files
   with gzip.open(files[0],'rt') as f:
    data=list(csv.DictReader(f,delimiter='\t'))
   metrics.append({'run_accession':a,'arm':arm,'final_pairs':len(data),'joint_exemplars':sum(x['seq_id']==x['prim_align_dup_exemplar']==x['sim_dup_exemplar'] for x in data)})
 (d/'metrics.json').write_text(json.dumps(metrics,indent=2))
 # Preserve task commands/logs before releasing only this batch's generated work files.
 with tarfile.open(d/'task-diagnostics.tar.gz','w:gz') as tar:
  for p in work.glob('*/*/.command*'):
   if p.is_file() and not p.is_symlink():tar.add(p,arcname=str(p.relative_to(work)),recursive=False)
 uri=f's3://{bucket}/runs/amino-acid-20261021/heldout-batches/{d.name}/'
 subprocess.run(['aws',*profile,'s3','cp',str(out)+'/',uri+'output/','--recursive','--only-show-errors'],check=True)
 for name in ['manifest.tsv','metrics.json','task-diagnostics.tar.gz','nextflow.log','run.log']:
  subprocess.run(['aws',*profile,'s3','cp',str(d/name),uri+name,'--only-show-errors'],check=True)
 result={'samples':ids,'batch':batch,'wall_seconds':time.monotonic()-start,'private_results':uri,'metrics_rows':len(metrics)}
 (d/'complete.json').write_text(json.dumps(result,indent=2))
 subprocess.run(['aws',*profile,'s3','cp',str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
 shutil.rmtree(work);shutil.rmtree(incoming)
 print(json.dumps({'batch':batch,'samples':len(ids),'completed_total':len(completed()),'wall_seconds':result['wall_seconds']}),flush=True)

if __name__=='__main__':
 keys=json.loads(aws('s3api','list-objects-v2','--bucket',bucket,'--prefix','runs/amino-acid-20261021/heldout-batches/','--query','Contents[?ends_with(Key, `/complete.json`)].Key','--output','json')) or []
 prior=set()
 for key in keys:prior.update(json.loads(aws('s3','cp',f's3://{bucket}/{key}','-','--only-show-errors'))['samples'])
 (s/'prior-completed.json').write_text(json.dumps(sorted(prior)))
 print(json.dumps({'prior_committed_samples':len(prior),'batch_prefix':os.environ.get('AA_BATCH_PREFIX','s')}),flush=True)
 while True:
  done=completed()
  if len(done)==len(rows):break
  if datetime.datetime.now(datetime.timezone.utc)>deadline:print('Scheduling deadline reached',flush=True);break
  available=sorted(ready()-done,key=order.get)
  if not available:
   if (s/'producer-finished.json').exists():break
   time.sleep(30);continue
  # First eight confirm full workflow wiring; later batches amortize orchestration.
  size=8 if not done else 25
  producer_done=(s/'producer-finished.json').exists()
  if len(available)<size and not producer_done:time.sleep(30);continue
  batch=1+len(list((s/'batches').glob('*/complete.json')))
  process_batch(available[:size],batch)
  if (r/'consumer-stop-after-batch').exists():break
 print(json.dumps({'completed':len(completed()),'planned':len(rows)}),flush=True)

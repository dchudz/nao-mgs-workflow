"""Resume only the identified full-cohort validation session after an engine finalization failure."""
import datetime,hashlib,json,re,shutil,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-cohort';session='98a3b25c-949c-4eac-adee-c0ef24cb3f78'
assert not (d/'complete.json').exists()
scope=json.loads((d/'scope.json').read_text());assert scope['planned']==1000 and len(scope['completed_samples'])==1000 and not scope['missing_samples']
original_log=(d/'nextflow.log').read_text();assert 'Session UUID: '+session in original_log and 'java.util.ConcurrentModificationException' in original_log
history=[line for line in (r/'.nextflow/history').read_text().splitlines() if session in line];assert len(history)==1 and 'validation-cohort/work' in history[0] and '\tERR\t' in history[0]
failed=d/'failed-engine-finalization-20261022T070344';failed.mkdir()
for path in [d/'nextflow.log',d/'run.log',d/'output/trace.tsv',r/'cohort-validation-scheduled.log']:
 if path.exists():shutil.copy2(path,failed/path.name)
sha={}
# Stream hashes on the worker's Python3.9.
for path in [d/'manifest.tsv',*sorted(d.glob('*-pooled.tsv.gz'))]:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 sha[path.name]=h.hexdigest()
assert len(sha)==6
cmd=[str(r/'nextflow'),'-log',str(d/'nextflow.log'),'run','experiments/amino_acids/validation-pooled.nf','-resume',session,'-c','experiments/amino_acids/validation-cohort.config','-work-dir',str(d/'work'),'--manifest',str(d/'manifest.tsv'),'--ref_dir',str(r/'reference'),'--output',str(d/'output')]
(d/'resume-protocol.json').write_text(json.dumps({'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'session':session,'command':cmd,'input_sha256':sha,'reason':'Nextflow ConcurrentModificationException during successful BBMERGE task finalization; no search parameter or input change.'},indent=2))
with (d/'run.log').open('w') as f:subprocess.run(cmd,cwd=r,stdout=f,stderr=f,check=True,timeout=28800)
assert len(list((d/'output/COPY_HITS').glob('*_validation_hits.tsv.gz')))==5
with (d/'analysis.log').open('w') as f:subprocess.run(['python3',str(r/'summarize-native-validation.py'),str(d),'--nodes',str(r/'reference/results/taxonomy-nodes.dmp')],stdout=f,stderr=f,check=True)
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-cohort/'
aws=['aws','--profile','pipeline-codex','s3','cp']
subprocess.run(aws+[str(d)+'/',uri,'--recursive','--exclude','work/*','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_samples':1000,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resumed_session':session,'summary_complete':True},indent=2))
subprocess.run(aws+[str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
print('COHORT_VALIDATION_COMPLETE',flush=True)

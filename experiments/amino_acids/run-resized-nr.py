"""Resume the fixed full query pool on non-interruptible AMD capacity after parity."""
import datetime,json,subprocess,time,hashlib,os,signal
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/combined-validation'
bucket='pipeline-lab-839993661825-us-east-1'
names=[x.split()[0] for x in (d/'combined-queries.fna').read_text().splitlines() if x.startswith('>')]
assert len(names)==817744 and len(names)==len(set(names))
assert hashlib.sha256((d/'combined-queries.fna').read_bytes()).hexdigest()=='6aa28bedbaf612c3a7c5a895b24a3be4f18d2c8264a439846427d943b289c597'
cmd=[str(r/'protein/tools/diamond'),'blastx','--db',str(r/'protein/current-nr/nr'),'--query',str(d/'combined-queries.fna'),'--out',str(d/'nr.tsv'),'--sensitive','--threads','24','--block-size','2','--index-chunks','4','--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--outfmt','6','qseqid','sseqid','bitscore','evalue','length','pident','staxids']
(d/'search-protocol.json').write_text(json.dumps({'command':cmd,'query_sha256':hashlib.sha256((d/'combined-queries.fna').read_bytes()).hexdigest(),'queries':len(names),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
stop=datetime.datetime.fromisoformat('2026-10-22T14:30:00+00:00')
remaining=(stop-datetime.datetime.now(datetime.timezone.utc)).total_seconds();assert remaining>600
with (d/'nr-search.log').open('w') as log:
 process=subprocess.Popen(['/usr/bin/time','-v']+cmd,stdout=log,stderr=log,start_new_session=True)
 try:
  code=process.wait(timeout=remaining)
  if code:raise subprocess.CalledProcessError(code,cmd)
 except subprocess.TimeoutExpired:
  os.killpg(process.pid,signal.SIGTERM)
  try:process.wait(timeout=15)
  except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
  (d/'incomplete.json').write_text(json.dumps({'reason':'Runtime/deadline cap reached; partial search is not a completed validation.','stopped_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'timeout_seconds':remaining},indent=2))
  for name in ['nr-search.log','search-protocol.json','incomplete.json']:
   subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/name),f's3://{bucket}/runs/amino-acid-20261021/protein-nr-combined-incomplete/'+name,'--only-show-errors'],check=True)
  raise

subprocess.run(['python3',str(r/'protein-nr-summary.py'),'--nodes',str(r/'protein/current-nr/nodes.dmp'),'--hits',str(d/'nr.tsv'),'--viral-calls',str(d/'pilot-viral-pairs.tsv'),str(d/'development-viral-pairs.tsv'),str(d/'full-viral-pairs.tsv'),str(d/'heldout-viral-pairs.tsv'),'--output',str(d/'validated-pairs.tsv')],check=True)
subprocess.run(['python3',str(r/'summarize-nr-controls.py')],check=True)
subprocess.run(['python3',str(r/'summarize-seed-nr.py')],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',f's3://{bucket}/runs/amino-acid-20261021/protein-nr-combined/', '--recursive','--exclude','complete.json','--only-show-errors'],check=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),f's3://{bucket}/runs/amino-acid-20261021/protein-nr-combined/complete.json','--only-show-errors'],check=True)
print('COMBINED_NR_COMPLETE',flush=True)

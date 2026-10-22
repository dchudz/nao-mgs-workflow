"""Separate fixed validation panel; no extrapolation to unsearched candidates."""
import datetime,json,subprocess,time,hashlib,os,signal,shutil
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/bounded-validation'
bucket='pipeline-lab-839993661825-us-east-1'
assert (d/'protocol.json').exists()
for name in ['native-controls.tsv','seed-selection.tsv']:
 shutil.copy2(r/'protein/combined-validation'/name,d/name)
cmd=[str(r/'protein/tools/diamond'),'blastx','--db',str(r/'protein/current-nr/nr'),'--query',str(d/'queries.fna'),'--out',str(d/'nr.tsv'),'--sensitive','--threads','16','--block-size','2','--index-chunks','4','--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--outfmt','6','qseqid','sseqid','bitscore','evalue','length','pident','staxids']
(d/'search-protocol.json').write_text(json.dumps({'command':cmd,'query_sha256':hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest(),'queries':sum(line.startswith('>') for line in (d/'queries.fna').read_text().splitlines()),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
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
   subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/name),f's3://{bucket}/runs/amino-acid-20261021/bounded-nr-incomplete/'+name,'--only-show-errors'],check=True)
  raise

subprocess.run(['python3',str(r/'protein-nr-summary.py'),'--nodes',str(r/'protein/current-nr/nodes.dmp'),'--hits',str(d/'nr.tsv'),'--viral-calls',str(d/'selected-viral-pairs.tsv'),'--output',str(d/'validated-pairs.tsv')],check=True)
for name in ['summarize-nr-controls.py','summarize-seed-nr.py']:
 subprocess.run(['python3',str(r/name),str(d)],check=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',f's3://{bucket}/runs/amino-acid-20261021/bounded-nr-results/','--recursive','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),f's3://{bucket}/runs/amino-acid-20261021/bounded-nr-results/complete.json','--only-show-errors'],check=True)
print('BOUNDED_NR_COMPLETE',flush=True)

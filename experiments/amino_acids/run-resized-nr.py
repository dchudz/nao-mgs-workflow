"""Resume the fixed full query pool on non-interruptible AMD capacity after parity."""
import datetime,json,subprocess,time,hashlib,os,signal
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/combined-validation'
bucket='pipeline-lab-839993661825-us-east-1'
assert int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemTotal:')))>=240*1024**2,'Expected larger-memory worker'
assert os.cpu_count()==32,'Expected32configuredcores'
assert int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('SwapTotal:')))>=63*1024**2,'Expected64GiBswapfallback'
names=[x.split()[0] for x in (d/'combined-queries.fna').read_text().splitlines() if x.startswith('>')]
assert len(names)==817744 and len(names)==len(set(names))
assert hashlib.sha256((d/'combined-queries.fna').read_bytes()).hexdigest()=='6aa28bedbaf612c3a7c5a895b24a3be4f18d2c8264a439846427d943b289c597'
cmd=[str(r/'protein/tools/diamond'),'blastx','--db',str(r/'protein/current-nr/nr'),'--query',str(d/'combined-queries.fna'),'--out',str(d/'nr.tsv'),'--sensitive','--threads','32','--block-size','2','--index-chunks','4','--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--outfmt','6','qseqid','sseqid','bitscore','evalue','length','pident','staxids']
(d/'search-protocol.json').write_text(json.dumps({'command':cmd,'query_sha256':hashlib.sha256((d/'combined-queries.fna').read_bytes()).hexdigest(),'queries':len(names),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'instance_type':'r7a.8xlarge','memory_GiB':256,'threads':32,'worker_hourly_cost_usd':2.4344+(1280*.08+3000*.005+875*.04)/730,'recovery_from':'nr24-oom-20261022T034211'},indent=2))
stop=datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00')
remaining=(stop-datetime.datetime.now(datetime.timezone.utc)).total_seconds();assert remaining>600
with (d/'nr-search.log').open('w') as log:
 process=subprocess.Popen(['/usr/bin/time','-v']+cmd,stdout=log,stderr=log,start_new_session=True)
 try:
  until=time.monotonic()+remaining
  while True:
   try:
    code=process.wait(timeout=min(30,max(0,until-time.monotonic())));break
   except subprocess.TimeoutExpired:
    info={line.split(':',1)[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemTotal:','MemAvailable:','SwapTotal:','SwapFree:'))}
    info['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();info['diamond_RSS_KiB']=None
    children=Path(f'/proc/{process.pid}/task/{process.pid}/children')
    if children.exists():
     for child in children.read_text().split():
      status=Path(f'/proc/{child}/status')
      try:
       for line in status.read_text().splitlines():
        if line.startswith('VmRSS:'):info['diamond_RSS_KiB']=int(line.split()[1])
      except FileNotFoundError:pass
    with (d/'memory-metrics.jsonl').open('a') as metrics:metrics.write(json.dumps(info)+'\n')
    if time.monotonic()>=until:raise
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

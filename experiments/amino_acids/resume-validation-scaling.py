"""Time the already prepared cap inputs while recording shared-worker concurrency."""
import datetime,hashlib,json,os,signal,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-scaling'
assert not (d/'results.json').exists(), 'Do not overwrite an existing timed run'
protocol=json.loads((d/'protocol.json').read_text())
protocol['deadline']='2026-10-22T04:00:00+00:00'
protocol['scheduling_note']='Timing resumed from unchanged prepared inputs before the long sensitive BLAST diagnostic finished. Explicit Docker CPU shares8192 match an eight-thread Nextflow task; other tasks share the worker, so timings are not isolated service latency. Repeated cap20 and true search CPU are retained.'
(d/'protocol.json').write_text(json.dumps(protocol,indent=2))
query_counts={}
for cap in protocol['caps']:
 p=d/f'cap{cap}';preparation=json.loads((p/'preparation.json').read_text())
 assert hashlib.sha256((p/'queries.fna').read_bytes()).hexdigest()==preparation['query_sha256']
 query_counts[cap]=preparation['query_count']
stop=datetime.datetime.fromisoformat('2026-10-22T04:00:00+00:00')
def run(cmd,log,cwd=None):
 remaining=(stop-datetime.datetime.now(datetime.timezone.utc)).total_seconds();assert remaining>60
 with log.open('w') as f:
  p=subprocess.Popen(cmd,cwd=cwd,stdout=f,stderr=f,start_new_session=True)
  try:
   code=p.wait(timeout=remaining)
   if code:raise subprocess.CalledProcessError(code,cmd)
  except subprocess.TimeoutExpired:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   raise
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
fmt='6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send'
results=[]
for index,cap in enumerate(protocol['execution_order']):
 b=d/f'run{index+1}-cap{cap}';b.mkdir(exist_ok=True)
 # GNU time is inside the container, measuring the search process rather than the Docker client.
 script='set -euo pipefail; /usr/bin/time -v -o "$1" blastn -db "$2" -perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads 8 -query "$3" -outfmt "$4" | gzip > "$5"'
 cmd=['docker','run','--rm','--name',f'aa-validation-scaling-{index+1}','--network','none','--cpu-shares','8192','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{r}:{r}',image,'bash','-c',script,'benchmark',str(b/'blast.time'),str(r/'scratch/blast_db_190b4ae1/blast_db'),str(d/f'cap{cap}/queries.fna'),fmt,str(b/'hits.tsv.gz')]
 start=time.monotonic()
 try:run(cmd,b/'run.log')
 finally:subprocess.run(['docker','stop','-t','10',f'aa-validation-scaling-{index+1}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 result={'cap':cap,'queries':query_counts[cap],'wall_seconds':time.monotonic()-start,'command':cmd,'run':index+1}
 (b/'resources.json').write_text(json.dumps(result,indent=2));results.append(result)
 (d/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='command'}),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-scaling/','--recursive','--exclude','*/work/*','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),'s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-scaling/complete.json','--only-show-errors'],check=True)

"""Compare BLAST query threading with its auto mode on the exact frozen cap20 input."""
import collections,datetime,gzip,hashlib,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');source=r/'cohort-storage/validation-scaling';d=r/'cohort-storage/blast-threading';d.mkdir(exist_ok=True)
assert not (d/'results.json').exists(),'Do not overwrite a prior trial'
query=source/'cap20/queries.fna';preparation=json.loads((source/'cap20/preparation.json').read_text());assert hashlib.sha256(query.read_bytes()).hexdigest()==preparation['query_sha256']
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'query_count':preparation['query_count'],'query_sha256':preparation['query_sha256'],'order':[1,0,1],'modes':{'0':'Pinned BLAST2.17 auto split by database or queries','1':'Explicit split by queries'},'same_inputs':'Unchanged cap20 query bytes, full core_nt, pinned image, scoring/search/output parameters and8threads.','resources':'Same primary worker; CPUshares8192,400GiB container memory ceiling on all new trials; two-hour cap per trial. Wait until the matched protein-NT diagnostic completes. No extra instance.','comparison':'Compare exact alignment multisets with completed cap20 run5 and between query-threaded repeats. No claim of faster identical behavior if alignment rows differ.','limits':'Shared worker and warm reference cache. Recorded wall time is not isolated latency or incremental billing.'}
(d/'protocol.json').write_text(json.dumps(protocol,indent=2));uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/blast-threading/'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'protocol.json'),uri+'protocol.json','--only-show-errors'],check=True)
while not (r/'cohort-storage/protein-nt-panel/complete.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime.fromisoformat('2026-10-23T06:00:00+00:00'),'Input wait exceeds budget reserve'
 time.sleep(30)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
fmt='6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send'
results=[]
for index,mode in enumerate(protocol['order'],1):
 b=d/f'run{index}-mode{mode}';b.mkdir();name=f'aa-blast-threading-{index}'
 script='set -euo pipefail; /usr/bin/time -v -o "$1" blastn -db "$2" -perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads 8 -mt_mode "$6" -query "$3" -outfmt "$4" | gzip > "$5"'
 cmd=['docker','run','--rm','--name',name,'--network','none','--cpu-shares','8192','--memory','400g','--memory-swap','400g','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{r}:{r}',image,'bash','-c',script,'benchmark',str(b/'blast.time'),str(r/'scratch/blast_db_190b4ae1/blast_db'),str(query),fmt,str(b/'hits.tsv.gz'),str(mode)]
 start=time.monotonic();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
 try:
  with (b/'run.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=7200)
 finally:subprocess.run(['docker','stop','-t','10',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 result=dict(run=index,mode=mode,queries=preparation['query_count'],started_utc=started,wall_seconds=time.monotonic()-start,command=cmd)
 (b/'resources.json').write_text(json.dumps(result,indent=2));results.append(result);(d/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='command'}),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri,'--recursive','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2));subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
print('BLAST_THREADING_COMPLETE',flush=True)

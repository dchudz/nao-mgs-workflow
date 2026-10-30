"""Matched nucleotide searches on a frozen sample-balanced protein-candidate panel."""
import datetime,hashlib,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/protein-nt-panel';stop=datetime.datetime.fromisoformat('2026-10-23T12:00:00+00:00')
assert hashlib.sha256((d/'queries.fna').read_bytes()).hexdigest()==json.loads((d/'query-provenance.json').read_text())['query_sha256']
while not (r/'cohort-storage/validation-scaling/complete.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Sampling benchmark did not complete within the budget safeguard'
 time.sleep(30)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
for task in ['megablast','blastn']:
 assert not (d/f'{task}.resources.json').exists(),'Refuse to replace completed measurements'
 name='aa-protein-nt-panel-'+task
 search=['blastn','-task',task,'-db',str(r/'scratch/blast_db_190b4ae1/blast_db'),'-query',str(d/'queries.fna'),'-out',str(d/f'{task}.tsv'),'-perc_identity','60','-max_hsps','5','-num_alignments','250','-qcov_hsp_perc','30','-num_threads','8','-outfmt','6 qseqid sseqid staxids bitscore evalue length pident qcovs']
 cmd=['docker','run','--rm','--name',name,'--network','none','--cpu-shares','8192','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{r}:{r}',image,'/usr/bin/time','-v','-o',str(d/f'{task}.time')]+search
 started=datetime.datetime.now(datetime.timezone.utc).isoformat();start=time.monotonic();remaining=(stop-datetime.datetime.now(datetime.timezone.utc)).total_seconds();assert remaining>600
 try:
  with (d/f'{task}.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=min(21600,remaining))
 finally:subprocess.run(['docker','stop','-t','10',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 row=dict(task=task,queries=464,pairs=232,started_utc=started,wall_seconds=time.monotonic()-start,command=cmd)
 (d/f'{task}.resources.json').write_text(json.dumps(row,indent=2));print(json.dumps({k:v for k,v in row.items() if k!='command'}),flush=True)
uri='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/protein-nt-panel/'
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/',uri,'--recursive','--exclude','complete.json','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps(dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),queries=464,pairs=232),indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),uri+'complete.json','--only-show-errors'],check=True)
print('PROTEIN_NT_PANEL_COMPLETE',flush=True)

"""Original query preparation and BLAST parameters, varying only the species cap."""
import csv,datetime,gzip,hashlib,json,os,signal,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-scaling';d.mkdir(exist_ok=True)
samples=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200]
found={a:p.parent for p in (r/'cohort-storage/batches').glob('*/complete.json') for a in json.loads(p.read_text())['samples']}
assert all(a in found for a in samples)
source=d/'baseline200.tsv.gz';header=None
with gzip.open(source,'wb') as out:
 for a in samples:
  files=list((found[a]/'output/marked').glob(f'*_{a}_nt_duplicate_reads.tsv.gz'));assert len(files)==1
  with gzip.open(files[0],'rb') as f:
   h=f.readline()
   if header is None:header=h;out.write(h)
   else:assert h==header
   for chunk in iter(lambda:f.read(1024*1024),b''):out.write(chunk)
(d/'manifest.tsv').write_text(f'id\tpath\nbaseline200\t{source}\n')
stop=datetime.datetime.fromisoformat('2026-10-22T02:00:00+00:00')
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
caps=[1,5,20,100]
protocol={'samples':samples,'caps':caps,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'input_scope':'First200 baseline NT24 sample tables, duplicate marking preserved within each sample. Species sampling is across the pooled200, not per sample.','query_preparation':'Original species partition, hash downsample with exemplar restriction, extract, merge/join, FASTA modules.','timed_scope':'Original BLAST search parameters against existing comprehensive core_nt plus gzip; excludes reference download and query preparation, which are separately recorded.','execution_order':[20,1,100,5,20],'concurrency':'Native cohort processing may share this worker; no claim of isolated service latency. Warm reference pages expected; repeated20 estimates temporal variation.','deadline':stop.isoformat(),'reference':'Existing297GiB core_nt database; no rebuild or bespoke reduction.'}
(d/'protocol.json').write_text(json.dumps(protocol,indent=2))
query_counts={};sets={}
for cap in caps:
 b=d/f'cap{cap}';b.mkdir(exist_ok=True)
 cmd=[str(r/'nextflow'),'-log',str(b/'nextflow.log'),'run','experiments/amino_acids/validation-scaling-queries.nf','-c','experiments/amino_acids/validation-scaling.config','-work-dir',str(b/'work'),'--manifest',str(d/'manifest.tsv'),'--ref_dir',str(r/'reference'),'--output',str(b/'output'),'--validation_n_sample',str(cap)]
 start=time.monotonic();run(cmd,b/'prepare.log',r)
 paths=list((b/'output/COPY_QUERIES').glob('*queries.fasta.gz'));assert len(paths)==1
 fasta=b/'queries.fna'
 with gzip.open(paths[0],'rb') as f:fasta.write_bytes(f.read())
 ids=[x.split()[0][1:] for x in fasta.read_text().splitlines() if x.startswith('>')];assert len(ids)==len(set(ids));sets[cap]=set(ids)
 query_counts[cap]=len(ids)
 (b/'preparation.json').write_text(json.dumps({'wall_seconds':time.monotonic()-start,'query_count':len(ids),'query_sha256':hashlib.sha256(fasta.read_bytes()).hexdigest(),'command':cmd},indent=2))
 print(json.dumps({'prepared_cap':cap,'queries':len(ids)}),flush=True)
assert all(sets[a]<=sets[b] for a,b in zip(caps,caps[1:])), 'Hash-selected queries must nest across caps'
while not (r/'cohort-storage/seed-diagnostic200/blastn.resources.json').exists():
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'Existing BLAST diagnostic did not finish before benchmark deadline'
 time.sleep(30)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
fmt='6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send'
results=[]
for index,cap in enumerate(protocol['execution_order']):
 b=d/f'run{index+1}-cap{cap}';b.mkdir(exist_ok=True)
 # GNU time is inside the container, measuring the search process rather than the Docker client.
 script='set -euo pipefail; /usr/bin/time -v -o "$1" blastn -db "$2" -perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads 8 -query "$3" -outfmt "$4" | gzip > "$5"'
 cmd=['docker','run','--rm','--name',f'aa-validation-scaling-{index+1}','--network','none','-v','/usr/bin/time:/usr/bin/time:ro','-v',f'{r}:{r}',image,'bash','-c',script,'benchmark',str(b/'blast.time'),str(r/'scratch/blast_db_190b4ae1/blast_db'),str(d/f'cap{cap}/queries.fna'),fmt,str(b/'hits.tsv.gz')]
 start=time.monotonic()
 try:run(cmd,b/'run.log')
 finally:subprocess.run(['docker','stop','-t','10',f'aa-validation-scaling-{index+1}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 result={'cap':cap,'queries':query_counts[cap],'wall_seconds':time.monotonic()-start,'command':cmd,'run':index+1}
 (b/'resources.json').write_text(json.dumps(result,indent=2));results.append(result)
 (d/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='command'}),flush=True)
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-scaling/','--recursive','--exclude','*/work/*','--only-show-errors'],check=True)
(d/'complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(d/'complete.json'),'s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/validation-scaling/complete.json','--only-show-errors'],check=True)

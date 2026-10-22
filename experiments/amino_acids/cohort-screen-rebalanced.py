"""Extract preselected public SRA blocks, verify pairs, and run matched screens."""
import concurrent.futures,subprocess,time,json,hashlib,sys,os,shutil,traceback
from pathlib import Path
r=Path('/home/ec2-user/experiment');sys.path.insert(0,str(r/'experiments/amino_acids'))
from union_pairs import pairs,union
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/rust-tools@sha256:e69b93d9a1a0ae41c59ec4c0a6a0841f3a7f4871fc2f9b333eaa094563612270'
def run(row):
 a=row['run_accession'];d=r/'cohort'/row['split']/a;d.mkdir(parents=True,exist_ok=True)
 if (d/'complete.json').exists():return json.loads((d/'complete.json').read_text())
 # Do not reuse partial FASTQ or outputs after failed extraction/screening.
 for name in ['reads.fastq','nt.fastq','nt20.fastq','nt22.fastq','aa12.fastq','union.fastq','rescued.fastq']:
  p=d/name
  if p.exists():raise RuntimeError(f'Partial output needs investigation: {p}')
 timings={};start=time.monotonic();raw=d/'reads.fastq'
 src=r/'sra'/a
 source=str(src) if src.exists() else f'https://sra-pub-run-odp.s3.amazonaws.com/sra/{a}/{a}'
 with raw.open('xb') as out:
  for i,(lo,hi) in enumerate(row['spot_ranges']):
   cmd=[str(r/'sratoolkit.3.4.1-ubuntu64/bin/fastq-dump'),'--split-spot','--skip-technical','--minReadLen','1','--defline-seq','@$ac.$si/$ri','--defline-qual','+','-N',str(lo),'-X',str(hi),'--stdout',source]
   part=d/f'block-{i}.partial'
   for attempt in range(3):
    try:
     with part.open('wb') as chunk,(d/f'extract-{i}-attempt-{attempt}.log').open('w') as log:subprocess.run(cmd,stdout=chunk,stderr=log,check=True,timeout=900)
     observed=0
     for key,data in pairs(part):
      assert key[0]==f'@{a}.{lo+observed}'.encode(),(key[0],lo+observed)
      observed+=1
     assert observed==hi-lo+1,(observed,hi-lo+1)
     break
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired):
     if attempt==2:raise
     time.sleep(attempt+1)
   with part.open('rb') as chunk:
    for data in iter(lambda:chunk.read(1024*1024),b''):out.write(data)
   part.unlink()
 timings['extract_seconds']=time.monotonic()-start
 h=hashlib.sha256();count=0
 for key,data in pairs(raw):h.update(data);count+=1
 expected=sum(hi-lo+1 for lo,hi in row['spot_ranges'])
 assert count==expected,(count,expected)
 commands={}
 for k in [24,22,20]:
  name='nt' if k==24 else f'nt{k}'
  idx=r/'reference/results/virus-genomes-masked.nucleaze.bin' if k==24 else r/f'reference/virus.nt{k}.bin'
  commands[name]=[str(r/'assets/nucleaze'),'--in',str(raw),'--interinput','--outm',str(d/f'{name}.fastq'),'--binref',str(idx),'--k',str(k),'--minhits','1','--canonical','--threads','2']
 commands['aa12']=[str(r/'rust-tools/target-speed/release/screen_amino_acids'),'screen','--index',str(r/'reference/virus.aa12.bin'),'--reads',str(raw),'--matched',str(d/'aa12.fastq'),'--threads','2']
 for name,cmd in commands.items():
  t=time.monotonic()
  with (d/f'{name}.log').open('w') as log:subprocess.run(['/usr/bin/time','-f','{\"wall_seconds\":%e,\"user_seconds\":%U,\"system_seconds\":%S,\"peak_rss_kib\":%M}','-o',str(d/f'{name}-resources.json')]+cmd,stdout=log,stderr=log,check=True,timeout=1800)
  timings[name+'_seconds']=time.monotonic()-t
 stats=union(d/'nt.fastq',d/'aa12.fastq',d/'union.fastq',d/'rescued.fastq')
 result={**row,**stats,'observed_pairs':count,'input_bytes':raw.stat().st_size,'input_sha256':h.hexdigest(),'timings':timings,'source':source,'execution_mode':os.environ.get('AA_SCREEN_EXECUTION','eight separate worker processes')}
 result['resources']={name:json.loads((d/f'{name}-resources.json').read_text()) for name in commands}
 (d/'complete.json').write_text(json.dumps(result,indent=2))
 if row['split']=='heldout':raw.unlink() # Retain selected pairs, hashes, manifest and public source; discard only generated scratch.
 print(json.dumps({'run':a,'pairs':count,**stats,'timings':timings}),flush=True)
 return result
def safe_run(row):
 a=row['run_accession'];d=r/'cohort'/row['split']/a
 try:
  assert shutil.disk_usage(r).free>25*1024**3,'Less than25GiB scratch space remains'
  result=run(row)
  if row['split']=='heldout':
   uri=f's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/cohort-screens/{a}/'
   subprocess.run(['aws','s3','cp',str(d)+'/',uri,'--recursive','--exclude','reads.fastq','--exclude','*.partial','--exclude','complete.json','--exclude','failure.json','--profile','pipeline-codex','--only-show-errors'],check=True)
   subprocess.run(['aws','s3','cp',str(d/'complete.json'),uri+'complete.json','--profile','pipeline-codex','--only-show-errors'],check=True)
  return result
 except Exception as exc:
  d.mkdir(parents=True,exist_ok=True)
  failure={'run_accession':a,'error':str(exc),'traceback':traceback.format_exc()}
  (d/'failure.json').write_text(json.dumps(failure,indent=2))
  for candidate in [d/'reads.fastq',*d.glob('*.partial')]:
   if candidate.exists():candidate.unlink() # Incomplete generated scratch, never a public source object.
  print(json.dumps({'run':a,'status':'failed','error':str(exc)}),flush=True)
  return failure

if __name__=='__main__':
 rows=json.loads(Path(sys.argv[1]).read_text())
 with concurrent.futures.ProcessPoolExecutor(int(sys.argv[2])) as pool:results=list(pool.map(safe_run,rows))
 Path(str(sys.argv[1])+'.results.json').write_text(json.dumps(results,indent=2))
 print(json.dumps({'completed':sum('error' not in x for x in results),'failed':sum('error' in x for x in results)}),flush=True)

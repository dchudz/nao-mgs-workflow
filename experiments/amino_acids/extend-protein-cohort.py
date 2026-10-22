"""Extend viral protein comparison to frozen samples 201–1000; no new NR search.

Uses existing AA12 FASTP files and unchanged contaminant/query/search algorithms.
Each additional 200-sample query pool remains separate. The original first200
pool, NR queries and results are never rewritten. Run --parity-only first.
"""
import argparse,csv,datetime,gzip,hashlib,json,os,shutil,subprocess,sys,tarfile,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');E=R/'experiments/amino_acids';S=R/'cohort-storage';D=S/'protein-extension800';D.mkdir(exist_ok=True)
sys.path.insert(0,str(E));from protein_queries import prepare
p=argparse.ArgumentParser();p.add_argument('--source-commit',required=True);p.add_argument('--parity-only',action='store_true');a=p.parse_args()
assert len(a.source_commit)==40 and all(c in '0123456789abcdef' for c in a.source_commit)
AWS=['aws','--profile','pipeline-codex','--region','us-east-1'];BUCKET='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
order=[x['run_accession'] for x in json.loads((R/'heldout-ranges.json').read_text())];assert len(order)==len(set(order))==1000
batches={s:m.parent for m in (S/'batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']};assert set(order)<=set(batches)
limit=time.monotonic()+14*3600
ABSOLUTE=datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00')
def now():return datetime.datetime.now(datetime.timezone.utc)
def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2)+'\n')
def run(cmd,log,cwd=R):
 remaining=min(limit-time.monotonic(),(ABSOLUTE-now()).total_seconds());assert remaining>0,'Budget execution window exhausted'
 with Path(log).open('w') as f:subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=remaining)
def backup(directory,prefix):
 skips=['work/*','*/work/*','upload.log','*/upload.log','verification.log','*/verification.log','backup-verification.json','*/backup-verification.json']
 excludes=[str(f.relative_to(directory)) for f in directory.rglob('*') if f.name in {'work','upload.log','verification.log','backup-verification.json'}]
 excludes+=['upload.log','verification.log','backup-verification.json']
 args=[x for pat in skips for x in ['--exclude',pat]]
 run(AWS+['s3','cp',str(directory)+'/',BUCKET+prefix+'/','--recursive',*args,'--only-show-errors'],directory/'upload.log')
 check=[x for item in sorted(set(excludes)) for x in ['--exclude-dir',item]]
 run(['python3',R/'verify-private-prefix.py',directory,'runs/amino-acid-20261021/'+prefix,*check,'--output',directory/'backup-verification.json'],directory/'verification.log')
def marker(directory,prefix,data):
 write(directory/'complete.json',dict(data,private_backups_verified=True))
 subprocess.run(AWS+['s3','cp',str(directory/'complete.json'),BUCKET+prefix+'/complete.json','--only-show-errors'],check=True)
def sources(ids):return [(sample,batches[sample]/'output/fastp'/f'{sample}_aa12_fastp.fastq.gz') for sample in ids]
def prepared_queries(ids,out):
 out.mkdir(exist_ok=True);q=out/'queries';q.mkdir(exist_ok=True)
 manifest=out/'manifest.tsv'
 with manifest.open('w') as f:
  f.write('id\tpath\n')
  for sample,path in sources(ids):assert path.is_file();f.write(f'{sample}_aa12\t{path}\n')
 write(out/'input-hashes.json',{sample:{'path':str(path),'sha256':digest(path),'bytes':path.stat().st_size} for sample,path in sources(ids)})
 run([R/'nextflow','-log',out/'nextflow.log','run','experiments/amino_acids/protein-contaminants-batch.nf','-c','experiments/amino_acids/protein-cohort.config','-work-dir',out/'work','--manifest',manifest,'--ref_dir',R/'reference','--output',out/'output'],out/'preparation.log')
 counts={}
 for sample in ids:
  source=out/'output/OTHER'/f'{sample}_aa12_other_bowtie2_unmapped.fastq.gz'
  counts[sample]=prepare(source,q/f'{sample}.fna',q/f'{sample}.tsv.gz',sample+'_')
 with tarfile.open(out/'task-diagnostics.tar.gz','w:gz') as tar:
  for f in (out/'work').glob('*/*/.command*'):
   if f.is_file() and not f.is_symlink():tar.add(f,arcname=str(f.relative_to(out/'work')))
 with (out/'all.fna').open('xb') as fasta,gzip.open(out/'all.tsv.gz','wt') as meta,gzip.open(out/'all.members.tsv.gz','wt') as members:
  for i,sample in enumerate(ids):
   with (q/f'{sample}.fna').open('rb') as f:shutil.copyfileobj(f,fasta)
   for suffix,outfile in [('.tsv.gz',meta),('.fna.members.tsv.gz',members)]:
    with gzip.open(q/(sample+suffix),'rt') as f:
     header=next(f)
     if i==0:outfile.write(header)
     shutil.copyfileobj(f,outfile)
 write(out/'preparation-complete.json',{'samples':ids,'counts':counts,'query_sha256':digest(out/'all.fna'),'completed_utc':now().isoformat()})
 return counts
# The single-sample characterization uses an original200 sample with a known
# read-header-sensitive contaminant case, preserving its exact original input.
parity=D/'preparation-parity'
if not (parity/'complete.json').exists():
 sample='SRR35904031';assert sample in order[:200];counts=prepared_queries([sample],parity)
 compared={}
 for suffix in ['.tsv.gz','.fna.members.tsv.gz']:
  old=S/'protein200/queries'/(sample+suffix)
  assert old.exists(),old
  with gzip.open(old,'rb') as f:left=f.read()
  with gzip.open(parity/'queries'/(sample+suffix),'rb') as f:right=f.read()
  assert left==right,('Preparation parity failed',suffix)
  compared[suffix]=hashlib.sha256(right).hexdigest()
 write(parity/'complete.json',{'sample':sample,'source_commit':a.source_commit,'completed_utc':now().isoformat(),'exact_decompressed_query_and_membership_parity':True,'files':compared,'counts':counts})
 print(json.dumps({'preparation_parity':True,'sample':sample}),flush=True)
if a.parity_only:raise SystemExit
assert json.loads((parity/'complete.json').read_text())['exact_decompressed_query_and_membership_parity']
protocol=D/'protocol.json';assert not protocol.exists(),'Extension already launched; inspect saved checkpoints before resuming'
identity=json.loads(subprocess.check_output(AWS+['sts','get-caller-identity'],text=True));assert identity['Account']=='839993661825'
source_files=['protein-contaminants-batch.nf','protein-cohort.config','local.config','protein_queries.py','union_pairs.py','protein-viral-summary.py','analyze-protein200.py','extend-protein-cohort.py']
exe=R/'protein/tools/diamond';assert digest(exe)=='c32782ed588593f73c41c57c0d013b7b8a5ac6d78e48d90a63985cbcf148d807'
write(protocol,{'source_commit':a.source_commit,'started_utc':now().isoformat(),'sample_selection':'Positions200:1000 of original frozen SHA-ranked1000sample manifest, four consecutive200sample query pools','samples':order[200:],'heldout_manifest_sha256':digest(R/'heldout-ranges.json'),'source_sha256':{n:digest(E/n) for n in source_files},'diamond_sha256':digest(exe),'viral_database_sha256':digest(R/'protein/viral-alignment.dmnd'),'database_bytes':(R/'protein/viral-alignment.dmnd').stat().st_size,'NR_scope':'No NR searches added. Originalfirst200NR queries/results unchanged. Supplemental calls unvalidated.','worker':'Existing private primary x2iedn.4xlarge; no new workers','budget':'Existing500USD cap and absolute worker guard unchanged; max14h driver budget and08UTC absolute cutoff. Additional occupied-primary planningreserve60USD fits within already-budgeted worker lifetime.','scientific_parameters':'Unchanged existing AA12FASTP bytes, human/other filtering, exact per-sample dedup, DIAMOND very-sensitive b2c4 minORF1 code1 E1e-5 qcover30 top10 min25aa and0.9near-bestLCA','comparability':'Prespecified separate200sample query pools. Do not assume cross-pool raw alignment equivalence. Original200 remains its own unchanged pool.'})
completed=[]
try:
 for i in range(4):
  ids=order[200+i*200:400+i*200];out=D/f'group{i+2:02d}';assert not out.exists(),out
  assert shutil.disk_usage(S).free>120*2**30,'Insufficient scratch reserve'
  prepared_queries(ids,out);write(out/'samples.json',ids)
  cmd=[exe,'blastx','--db',R/'protein/viral-alignment','--query',out/'all.fna','--out',out/'viral.tsv','--very-sensitive','--threads','16','--block-size','2','--index-chunks','4','--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--outfmt','6','qseqid','sseqid','bitscore','evalue','length','pident','qlen','qstart','qend','sstart','send','qframe']
  write(out/'search-protocol.json',{'source_commit':a.source_commit,'command':list(map(str,cmd)),'query_sha256':digest(out/'all.fna'),'started_utc':now().isoformat(),'samples':ids,'query_mates':2*sum(v['exact_unique_pairs'] for v in json.loads((out/'preparation-complete.json').read_text())['counts'].values())})
  run(['/usr/bin/time','-v',*cmd],out/'viral-search.log')
  run(['python3',E/'protein-viral-summary.py',out/'viral.tsv',out/'all.tsv.gz',out/'viral-pairs.tsv'],out/'viral-summary.log')
  for arm in ['nt','nt20']:
   run(['python3',E/'analyze-protein200.py','--protein',out,'--output',out/f'comparison-{arm}','--samples-json',out/'samples.json','--native-arm',arm],out/f'comparison-{arm}.log')
  write(out/'search-complete.json',{'samples':ids,'source_commit':a.source_commit,'completed_utc':now().isoformat(),'NR_run':False,'broad_reference_status':'unvalidated'})
  backup(out,'protein-extension800/'+out.name)
  marker(out,'protein-extension800/'+out.name,json.loads((out/'search-complete.json').read_text()))
  shutil.rmtree(out/'work')
  completed.extend(ids);write(D/'progress.json',{'completed_samples':len(completed),'samples':completed,'updated_utc':now().isoformat()});print(json.dumps({'supplemental_protein_samples':len(completed),'total_with_original200':len(completed)+200}),flush=True)
 backup(D,'protein-extension800')
 marker(D,'protein-extension800',{'source_commit':a.source_commit,'completed_utc':now().isoformat(),'samples':completed,'completed_samples':len(completed),'original200_unchanged':True,'NR_run':False})
except BaseException as exc:
 write(D/'incomplete.json',{'error':repr(exc),'completed_samples':len(completed),'recorded_utc':now().isoformat()});raise

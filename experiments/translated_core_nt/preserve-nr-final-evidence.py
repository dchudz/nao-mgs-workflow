"""Preserve all remaining NR experiment evidence before any worker retirement."""
import datetime,hashlib,io,json,os,subprocess,tarfile,tempfile
from pathlib import Path
r=Path('/home/ec2-user/experiment');protein=r/'protein';out=r/'nr-closure-final';out.mkdir(exist_ok=True)
assert (protein/'seed-nr-standalone/complete.json').exists() and (protein/'protein232-nr-evidence/complete.json').exists()
assert not subprocess.run(['pgrep','-x','diamond'],capture_output=True).stdout.strip(),'A measured DIAMOND search is still active'
archive=json.loads((protein/'reference-archive/complete.json').read_text());assert archive['verified'] and len(archive['files'])==1619
raw=json.loads((protein/'nr-query-call-full/backup-download-verification.json').read_text());assert raw['all_download_hashes_verified']
parity=json.loads((protein/'full-aggregation-parity.json').read_text());assert parity['byte_identical'] and parity['pair_rows']==406173
a=['aws','--profile','pipeline-codex','--region','us-east-1'];bucket='pipeline-lab-839993661825-us-east-1';prefix='runs/amino-acid-20261021/nr-closure-final/'
assert subprocess.check_output(a+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
skip={protein/'combined-validation/nr.tsv',protein/'nr-query-call-full/query-calls.tsv',protein/'nr-query-call-full/native-control-alignments.tsv',protein/'nr-query-call-full/seed-alignments.tsv'}
# Those four payloads have prior independent S3 download-SHA verification; retain their receipts in this archive.
paths=[];links=[]
for directory,names,files in os.walk(protein,followlinks=False):
 if Path(directory)==protein and 'current-nr' in names:names.remove('current-nr')
 for name in list(names):
  p=Path(directory)/name
  if p.is_symlink():links.append({'path':str(p.relative_to(r)),'target':os.readlink(p)});names.remove(name)
 for name in files:
  p=Path(directory)/name
  if p.is_symlink():links.append({'path':str(p.relative_to(r)),'target':os.readlink(p)})
  elif p not in skip:paths.append(p)
for name in ['nodes.dmp','names.dmp','nr.pal','nr.pjs','nr-prot-metadata.json','source-manifest.json']:
 p=protein/'current-nr'/name;assert p.is_file() and not p.is_symlink();paths.append(p)
for p in r.iterdir():
 if p.is_file() and not p.is_symlink() and p.suffix in ['.py','.cpp','.rs','.sh','.json','.jsonl','.tsv','.txt','.log','.nf']:paths.append(p)
paths=sorted(set(paths));manifest=[];target=out/'analysis-and-source.tar.gz';assert not target.exists(),'Preserve previous closure attempts'
class Reader:
 def __init__(self,f,h):self.f=f;self.h=h;self.bytes=0
 def read(self,n):
  b=self.f.read(n);self.h.update(b);self.bytes+=len(b);return b
with tarfile.open(target,'w:gz',compresslevel=1) as tar:
 for path in paths:
  stat=path.stat();h=hashlib.sha256();info=tar.gettarinfo(str(path),arcname=str(path.relative_to(r)))
  with path.open('rb') as f:
   reader=Reader(f,h);tar.addfile(info,reader)
  assert reader.bytes==stat.st_size
  mutable_root_log=path.parent==r and path.suffix in ['.log','.jsonl']
  if not mutable_root_log:assert path.stat().st_size==stat.st_size and path.stat().st_mtime_ns==stat.st_mtime_ns
  manifest.append({'path':str(path.relative_to(r)),'bytes':reader.bytes,'sha256':h.hexdigest(),'root_log_snapshot':mutable_root_log})
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
digest=sha(target);uri='s3://'+bucket+'/'+prefix+target.name;subprocess.run(a+['s3','cp',str(target),uri,'--only-show-errors'],check=True)
with tempfile.TemporaryDirectory(dir=out) as tmp:
 copy=Path(tmp)/target.name;subprocess.run(a+['s3','cp',uri,str(copy),'--only-show-errors'],check=True);assert copy.stat().st_size==target.stat().st_size and sha(copy)==digest
record={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'account':'839993661825','archive_uri':uri,'archive_bytes':target.stat().st_size,'archive_sha256':digest,'archive_download_hash_verified':True,'files':manifest,'excluded_symlinks':links,'previously_download_hash_verified_payloads':raw['files'],'frozen_NR_reference_archive_prefix':archive['prefix'],'frozen_NR_reference_objects':len(archive['files']),'scope':'All remaining regular NR scientific outputs, failed-attempt evidence, taxonomy files and root source/log snapshots. Four large already-verified raw/cache payloads remain in their original private S3 objects. Public reference volumes are in the separately verified Deep Archive snapshot. No worker state is changed.'};(out/'complete.json').write_text(json.dumps(record,indent=2)+'\n')
subprocess.run(a+['s3','cp',str(out/'complete.json'),'s3://'+bucket+'/'+prefix+'complete.json','--only-show-errors'],check=True);data=subprocess.check_output(a+['s3','cp','s3://'+bucket+'/'+prefix+'complete.json','-','--only-show-errors']);assert hashlib.sha256(data).hexdigest()==sha(out/'complete.json')
print(json.dumps({k:v for k,v in record.items() if k not in ['files','excluded_symlinks','previously_download_hash_verified_payloads']}),flush=True)

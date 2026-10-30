"""Prepare an immutable private evidence archive; upload/verify separately."""
import hashlib,json,tarfile
from pathlib import Path
r=Path(__file__).resolve().parent;a=r/'anchor-analysis';out=r/'anchor-result-20261023.tar'
files=[p for p in (a/'result-v1').rglob('*') if p.is_file()]+list(a.glob('*.json'))+[a/'anchor-query-calls.tsv.gz',r/'anchor-nr-result.html']
files += [r/n for n in ['analyze-anchor.py','restore-anchor-inputs.py','verify-anchor-result.py','write-anchor-report.py','preserve-anchor-result.py']]
binary=r.parent/'anchor-tools/source/target/release/mark_duplicates_similarity';files.append(binary)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
manifest={}
with tarfile.open(out,'w') as tar:
 for p in sorted(files):
  name=str(p.relative_to(r)) if p.is_relative_to(r) else 'tools/'+p.name
  assert name not in manifest;manifest[name]={'sha256':sha(p),'bytes':p.stat().st_size};tar.add(p,arcname=name)
receipt={'uri':'s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/anchor-analysis-20261023/'+out.name,'sha256':sha(out),'bytes':out.stat().st_size,'files':manifest,'original_inputs':'Original full and NR archive locations and checksums are recorded in the enclosed restoration receipts. No new sequence searches were run.'}
(a/'preservation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))

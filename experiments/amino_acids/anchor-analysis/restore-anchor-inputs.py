"""Restore only named, checksum-verified files from private experiment archives."""
import argparse,hashlib,json,tarfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('scope',choices=['nr','full']);a=p.parse_args();r=Path(__file__).resolve().parent;out=r/'anchor-analysis';out.mkdir(exist_ok=True)
if a.scope=='nr':
 receipt=json.loads((r/'nr-final-evidence-complete.json').read_text());archive=r/'anchor-nr-closure.tar.gz';digest=receipt['archive_sha256'];uri=receipt['archive_uri']
 wanted={'protein/current-nr/nodes.dmp','protein/current-nr/names.dmp','protein/compacted-validation-v2/validated-pairs.tsv','protein/compacted-validation-v2/competitor-clades.tsv','protein/compacted-validation-v2/query-call-cache.json','protein/combined-validation/search-protocol.json','protein/combined-validation/nr-search.log'}
else:
 receipt=json.loads((r/'primary-final-evidence/full.json').read_text());archive=r/'anchor-full.tar.gz';digest=receipt['archive']['sha256'];uri=receipt['archive']['uri']
 wanted={'full/protein/queries.tsv.gz','full/protein/queries.fna.members.tsv.gz','full/protein/viral-pairs.tsv','full/protein/native-comparison.tsv','full/protein/contaminants/OTHER/FULL_SRR39694979_other_bowtie2_unmapped.fastq.gz'}
 wanted.update(x['path'] for x in receipt['files'] if ('/marked/' in x['path'] and ('_nt_' in x['path'] or '_nt20_' in x['path'])) or x['path'].endswith(('_nt_fastp.fastq.gz','_nt20_fastp.fastq.gz')))
 wanted.update(x['path'] for x in receipt['files'] if '/RENAME_VIRUS_HITS/' in x['path'] or x['path']=='full/control/postscreen/FASTP/SRR39694979_fastp.fastq.gz')
 wanted.update(x['path'] for x in receipt['files'] if x['path'].endswith('.json') and ('screen' in x['path'] or 'protocol' in x['path']))
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
assert sha(archive)==digest,'Archive digest mismatch'
manifest={x['path']:x for x in receipt['files']};found=[]
with tarfile.open(archive,'r|gz') as tar:
 for member in tar:
  name=member.name.removeprefix('./')
  if name not in wanted:continue
  assert member.isfile();target=out/name;target.parent.mkdir(parents=True,exist_ok=True)
  h=hashlib.sha256()
  with tar.extractfile(member) as src,target.open('xb') as dst:
   for b in iter(lambda:src.read(8*1024**2),b''):h.update(b);dst.write(b)
  assert h.hexdigest()==manifest[name]['sha256'],name
  found.append(name)
assert set(found)==wanted,(wanted-set(found))
(out/(a.scope+'-restore-receipt.json')).write_text(json.dumps({'archive_uri':uri,'archive_sha256':digest,'verified_files':{n:manifest[n] for n in found}},indent=2)+'\n');print(json.dumps({'scope':a.scope,'verified_files':len(found)}))

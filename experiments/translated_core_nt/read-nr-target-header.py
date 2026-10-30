"""Preserve one exact native BLAST defline using previously verified PIN offsets."""
import hashlib,json,struct,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'protein/seed-nr-forced-double-index-trace/output-audit-v2'
index=json.loads((d/'human-target-location.json').read_text());target=index['targets'][0]
v=next(x for x in index['volumes'] if x['database']==target['database'])
oid=target['local_ordinal'];assert 0<=oid<v['sequences']
with Path(v['database']+'.pin').open('rb') as f:
    assert hashlib.sha256(f.read(v['header_bytes'])).hexdigest()==v['header_sha256']
    f.seek(v['header_bytes']+4*oid);start,end=struct.unpack('>II',f.read(8))
assert 0<end-start<1024**2
with Path(v['database']+'.phr').open('rb') as f:f.seek(start);data=f.read(end-start)
assert len(data)==end-start
output=d/'human-target-815174796.header.asn1';assert not output.exists();output.write_bytes(data)
cmd=['openssl','asn1parse','-inform','DER','-in',str(output)]
decoded=subprocess.run(cmd,check=True,capture_output=True,text=True)
(d/'human-target-815174796.header.txt').write_text(decoded.stdout)
receipt=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),global_ordinal=target['global_ordinal'],volume=v['database'],local_ordinal=oid,header_offset=start,header_bytes=len(data),header_sha256=hashlib.sha256(data).hexdigest(),decode_command=cmd,limits=['Exact stored defline and ASN.1 decoding, not inferred protein annotation.','PIN field layout follows the already-pinned DIAMOND 2.1.17 implementation. This retrieves one header, not a new homology search.'])
(d/'human-target-header-retrieval.json').write_text(json.dumps(receipt,indent=2)+'\n');print(decoded.stdout)

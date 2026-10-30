"""Compare generic gzip settings on one existing verified artifact, locally."""
import datetime as dt
import gzip
import hashlib
import json
import os
import platform
import statistics
import time
import zlib
from pathlib import Path

root=Path(__file__).resolve().parent;out=root/'compression-control';out.mkdir(exist_ok=True)
source=root/'interrupted-worker-completed-evidence.tar.gz';proof=json.loads((root/'interrupted-worker-completed-backup.json').read_text())
assert source.stat().st_size==proof['bytes'] and hashlib.file_digest(source.open('rb'),'sha256').hexdigest()==proof['sha256']
raw=out/'payload.tar';h=hashlib.sha256();size=0
with gzip.open(source,'rb') as inp,raw.open('wb') as dst:
    for block in iter(lambda:inp.read(8*1024**2),b''):dst.write(block);h.update(block);size+=len(block)
    dst.flush();os.fsync(dst.fileno())
expected=h.hexdigest();rows=[]
for trial,level in enumerate([9,1,6,6,1,9],1):
    target=out/f'level-{level}.tar.gz';start=time.perf_counter()
    with raw.open('rb') as inp,target.open('wb') as file:
        with gzip.GzipFile(filename='',mode='wb',fileobj=file,compresslevel=level,mtime=0) as dst:
            for block in iter(lambda:inp.read(8*1024**2),b''):dst.write(block)
        file.flush();os.fsync(file.fileno())
    seconds=time.perf_counter()-start;restored=hashlib.sha256();restored_size=0
    with gzip.open(target,'rb') as inp:
        for block in iter(lambda:inp.read(8*1024**2),b''):restored.update(block);restored_size+=len(block)
    assert restored.hexdigest()==expected and restored_size==size
    row={'trial':trial,'gzip_level':level,'compression_and_fsync_seconds':seconds,'compressed_bytes':target.stat().st_size,'compressed_sha256':hashlib.file_digest(target.open('rb'),'sha256').hexdigest(),'restored_payload_sha256':restored.hexdigest(),'verified':True}
    rows.append(row);(out/'partial-results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
summary=[]
for level in [1,6,9]:
    r=[x for x in rows if x['gzip_level']==level];assert len({x['compressed_sha256'] for x in r})==1
    summary.append({'gzip_level':level,'median_seconds':statistics.median(x['compression_and_fsync_seconds'] for x in r),'compressed_bytes':r[0]['compressed_bytes'],'restored_payload_verified':True})
result={'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source_archive':proof,'payload_bytes':size,'payload_sha256':expected,'python':platform.python_version(),'platform':platform.platform(),'machine':platform.machine(),'zlib_runtime':zlib.ZLIB_RUNTIME_VERSION,'trials':rows,'summary':summary,'limits':['One already completed artifact on the local computer, not an AWS or end-to-end pipeline benchmark.','Two observations per level in reversed order; no inferential timing claim.','Compression timing includes file fsync but excludes decompression verification.','No sequences were interpreted, no pipeline setting changed, and no existing archive was replaced.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'summary':summary,'payload_bytes':size}),flush=True)

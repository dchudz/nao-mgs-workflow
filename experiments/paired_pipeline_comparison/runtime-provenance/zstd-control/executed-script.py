"""Test an alternative lossless archive codec on the same fixed local payload."""
import compression.zstd as zstd
import datetime as dt
import hashlib
import json
import os
import platform
import time
from pathlib import Path

root=Path(__file__).resolve().parent;previous=json.loads((root/'compression-control/result.json').read_text());raw=root/'compression-control/payload.tar';out=root/'zstd-control';out.mkdir(exist_ok=True)
assert raw.stat().st_size==previous['payload_bytes'] and hashlib.file_digest(raw.open('rb'),'sha256').hexdigest()==previous['payload_sha256']
rows=[]
for trial in [1,2]:
    target=out/'payload.tar.zst';start=time.perf_counter()
    with raw.open('rb') as inp,target.open('wb') as file:
        with zstd.ZstdFile(file,mode='w',level=3) as dst:
            for block in iter(lambda:inp.read(8*1024**2),b''):dst.write(block)
        file.flush();os.fsync(file.fileno())
    elapsed=time.perf_counter()-start;h=hashlib.sha256();size=0
    with zstd.ZstdFile(target,'r') as inp:
        for block in iter(lambda:inp.read(8*1024**2),b''):h.update(block);size+=len(block)
    assert h.hexdigest()==previous['payload_sha256'] and size==previous['payload_bytes']
    r={'trial':trial,'codec':'zstandard','level':3,'compression_and_fsync_seconds':elapsed,'compressed_bytes':target.stat().st_size,'compressed_sha256':hashlib.file_digest(target.open('rb'),'sha256').hexdigest(),'restored_payload_sha256':h.hexdigest(),'verified':True};rows.append(r);print(json.dumps(r),flush=True)
assert len({r['compressed_sha256'] for r in rows})==1
result={'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source_archive':previous['source_archive'],'payload_bytes':previous['payload_bytes'],'payload_sha256':previous['payload_sha256'],'python':platform.python_version(),'platform':platform.platform(),'zstd_version':zstd.zstd_version,'trials':rows,'executed_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limits':['Same single local payload used in the gzip control; this is not an AWS or end-to-end pipeline benchmark.','Zstandard level 3 is one tested configuration, not an optimized or matched compression-ratio comparison.','Timing includes fsync and excludes subsequent decompression verification.','A Zstandard reader is required; deployment compatibility was not tested.','No original artifact, scientific input, or running pipeline was changed.']}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')

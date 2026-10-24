"""Verify a retained index's bytes, using read-only concurrent I/O first."""
import concurrent.futures
import datetime
import hashlib
import json
import os
from pathlib import Path
import sys
import time

path=Path('/home/ec2-user/experiment/protein/viral-alignment.dmnd')
expected='d154a0b68879119bbfd933edfcb8fe6a2d9cfffc7a1f28eb5535ce8d1bbaab04'
start=time.monotonic();before=path.stat();fd=os.open(path,os.O_RDONLY)
record={'worker':sys.argv[1],'path':str(path),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'bytes':before.st_size,'expected_sha256':expected,'mode':'Read-only storage initialization and checksum; no index rebuild or search-setting change'}
def read_region(lo):
    hi=min(lo+64*1024**2,before.st_size);pos=lo
    while pos<hi:
        block=os.pread(fd,min(1024**2,hi-pos),pos)
        if not block:raise RuntimeError('Unexpected end of file')
        pos+=len(block)
try:
    with concurrent.futures.ThreadPoolExecutor(32) as pool:
        list(pool.map(read_region,range(0,before.st_size,64*1024**2)))
finally:os.close(fd)
h=hashlib.sha256()
with path.open('rb') as f:
    for block in iter(lambda:f.read(16*1024**2),b''):h.update(block)
after=path.stat();record.update(sha256=h.hexdigest(),seconds=time.monotonic()-start,
    finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),passed=h.hexdigest()==expected and before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns)
Path('/home/ec2-user/paired-study-20261024/alignment-index-storage-check.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
assert record['passed']

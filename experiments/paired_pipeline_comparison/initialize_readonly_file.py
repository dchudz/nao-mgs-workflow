"""Read each byte of one immutable file with bounded concurrent I/O."""
import concurrent.futures
import datetime
import json
import os
from pathlib import Path
import sys
import time

path=Path('/home/ec2-user/experiment/scratch/kraken_db_9f172f33/hash.k2d')
out=Path('/home/ec2-user/paired-study-20261024/reference-read-initialization.json')
before=path.stat();size=before.st_size;start=time.monotonic()
record={'worker':sys.argv[1],'path':str(path),'bytes':size,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'mode':'read-only','concurrent_readers':32,'block_bytes':1024**2,'region_bytes':1024**3,
        'scope':'Storage initialization only. No alignment, classification, sequence processing, reference edits or search-setting changes.',
        'timing_limit':'Concurrent I/O intervention; first execution is not a steady-state timing measurement.'}
out.with_suffix('.started.json').write_text(json.dumps(record,indent=2)+'\n')
fd=os.open(path,os.O_RDONLY)
def read_region(bounds):
    lo,hi=bounds;pos=lo
    while pos<hi:
        block=os.pread(fd,min(1024**2,hi-pos),pos)
        if not block:raise RuntimeError('Unexpected end of immutable file')
        pos+=len(block)
    return hi-lo
try:
    regions=[(lo,min(lo+1024**3,size)) for lo in range(0,size,1024**3)]
    with concurrent.futures.ThreadPoolExecutor(32) as pool:
        total=0
        for n,count in enumerate(pool.map(read_region,regions),1):
            total+=count
            if n%16==0:print(json.dumps({'worker':sys.argv[1],'completed_regions':n,'total_regions':len(regions)}),flush=True)
    after=path.stat()
    assert total==size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
finally:os.close(fd)
record.update(completed=True,read_bytes=total,seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)

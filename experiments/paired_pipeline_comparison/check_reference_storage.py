"""Read-only checksum of the shared profiling index; record cache intervention."""
import datetime
import hashlib
import json
from pathlib import Path
import sys
import time

path=Path('/home/ec2-user/experiment/scratch/kraken_db_9f172f33/hash.k2d')
out=Path('/home/ec2-user/paired-study-20261024/reference-storage-check.json')
start=time.monotonic();before=path.stat();h=hashlib.sha256();size=0
record={'worker':sys.argv[1],'path':str(path),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'bytes':before.st_size,'purpose':'Read-only byte-identity verification and file-cache population on retained/restored storage',
        'timing_limit':'Concurrent index reading affects first-run wall times; include this storage work in experiment cost and do not call these steady-state deployment timings.'}
out.with_suffix('.started.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'worker':sys.argv[1],'started':True,'bytes':before.st_size}),flush=True)
with path.open('rb') as f:
    for block in iter(lambda:f.read(16*1024**2),b''):
        h.update(block);size+=len(block)
after=path.stat()
assert size==before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
record.update(sha256=h.hexdigest(),seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
out.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)

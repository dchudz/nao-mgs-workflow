"""Read and checksum an already-downloaded archive without modifying it."""
import datetime
import hashlib
import json
import sys
import time
from pathlib import Path

worker, sample = sys.argv[1:]
root = Path('/home/ec2-user/paired-study-20261024')
source = root/'inputs'/sample/sample
assert source.is_file() and source.name.startswith('SRR')
before = source.stat()
start = datetime.datetime.now(datetime.timezone.utc)
t = time.monotonic(); h = hashlib.sha256(); count = 0
with source.open('rb') as f:
    for block in iter(lambda:f.read(16*1024**2), b''):
        h.update(block); count += len(block)
after = source.stat()
assert count == before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
result = dict(worker=worker, sample=sample, path=str(source), bytes=count, sha256=h.hexdigest(),
              started_utc=start.isoformat(), finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              seconds=time.monotonic()-t, read_only=True,
              reason='Archive validator observed waiting on file pages; a sequential checksum verifies bytes and changes cache state without changing data or validation settings.',
              interpretation='Concurrent storage intervention. Its duration is not isolated pipeline time or a controlled speed measurement.')
target=root/(worker+'-'+sample+'-source-archive-storage-check.json')
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))

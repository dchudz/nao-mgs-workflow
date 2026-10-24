"""Hash immutable cache files with bounded read-only I/O; no database queries."""
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path
import sys
import time

root=Path('/home/ec2-user/experiment/scratch/blast_db_190b4ae1')
out=Path('/home/ec2-user/paired-study-20261024/validation-cache-storage-check.json')
files=sorted((p for p in root.iterdir() if p.is_file()),key=lambda p:p.stat().st_size,reverse=True)
start=time.monotonic()
record={'worker':sys.argv[1],'directory':str(root),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'files_expected':len(files),'bytes_expected':sum(p.stat().st_size for p in files),'concurrent_readers':32,
        'scope':'Read-only opaque-file checksums and storage initialization, not a validation search. No reference or method changes.',
        'timing_limit':'Concurrent I/O affects early stage times and remains part of experiment spending.'}
out.with_suffix('.started.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
def check(p):
    before=p.stat();h=hashlib.sha256();size=0
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block);size+=len(block)
    after=p.stat();assert size==before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
    return {'name':p.name,'bytes':size,'sha256':h.hexdigest()}
with concurrent.futures.ThreadPoolExecutor(32) as pool:
    rows=[]
    for future in concurrent.futures.as_completed([pool.submit(check,p) for p in files]):
        rows.append(future.result())
        if len(rows)%100==0:print(json.dumps({'worker':sys.argv[1],'files_checked':len(rows),'bytes_checked':sum(x['bytes'] for x in rows)}),flush=True)
assert {p.name for p in root.iterdir() if p.is_file()}=={r['name'] for r in rows}
rows.sort(key=lambda r:r['name']);manifest_hash=hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest()
record.update(completed=True,files=rows,manifest_sha256=manifest_hash,seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
out.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='files'}),flush=True)

"""Reclaim duplicate local core_nt storage without changing a byte or lookup path.

Only the immutable experiment reference and its production download_db cache are
eligible. Use the production cache lock, verify full byte equality, then replace
the cached file atomically with a hard link to the original reference file.
"""
import concurrent.futures,datetime,fcntl,hashlib,json,os,shutil,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');source=R/'reference/results/blast_db'
cache=R/'scratch'/('blast_db_'+hashlib.sha256(str(source).encode()).hexdigest()[:8])
assert cache.name=='blast_db_190b4ae1'
assert json.loads((R/'setup/reference-verification.json').read_text())['verified']
assert source.stat().st_dev==cache.stat().st_dev
source_files={str(p.relative_to(source)):p for p in source.rglob('*') if p.is_file()}
cache_files={str(p.relative_to(cache)):p for p in cache.rglob('*') if p.is_file()}
assert source_files.keys()==cache_files.keys()
start=time.monotonic();free=shutil.disk_usage(R).free
journal=R/'setup/reference-cache-links.jsonl'
def share(rel):
    src=source_files[rel];dst=cache_files[rel]
    assert not src.is_symlink() and not dst.is_symlink()
    before_src=src.stat();before_dst=dst.stat();assert before_src.st_size==before_dst.st_size
    h=hashlib.sha256()
    with src.open('rb') as a,dst.open('rb') as b:
        while True:
            x=a.read(8*1024**2);y=b.read(8*1024**2)
            assert x==y,('Reference cache differs',rel)
            if not x:break
            h.update(x)
    assert src.stat().st_mtime_ns==before_src.st_mtime_ns and dst.stat().st_mtime_ns==before_dst.st_mtime_ns
    already=before_src.st_ino==before_dst.st_ino
    if not already:
        temporary=dst.with_name(dst.name+'.verified-link')
        assert not temporary.exists()
        os.link(src,temporary);os.replace(temporary,dst)
    assert os.path.samefile(src,dst)
    return {'path':rel,'bytes':before_src.st_size,'sha256':h.hexdigest(),'already_shared':already}
with (cache.parent/(cache.name+'.lock')).open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    with journal.open('a') as log,concurrent.futures.ThreadPoolExecutor(4) as pool:
        rows=[]
        for row in pool.map(share,sorted(source_files)):
            rows.append(row);log.write(json.dumps(row)+'\n');log.flush()
receipt={'verified':True,'source':str(source),'cache':str(cache),'files':rows,
 'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seconds':time.monotonic()-start,
 'bytes_shared':sum(x['bytes'] for x in rows if not x['already_shared']),
 'observed_free_space_delta_bytes':shutil.disk_usage(R).free-free,
 'method':'Full byte comparison and SHA256 before atomic hard-link replacement, under the unchanged production cache lock. Reference bytes and all lookup paths preserved.'}
(R/'setup/reference-cache-links.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))

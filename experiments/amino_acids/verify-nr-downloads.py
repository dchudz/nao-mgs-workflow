"""Verify private NR backups by streaming downloads against recorded source SHA256."""
import datetime
import hashlib
import json
import subprocess
import time
from pathlib import Path

r = Path('/home/ec2-user/experiment')
d = r / 'protein/combined-validation'
cache = r / 'protein/nr-query-call-full'
m = json.loads((d / 'query-call-cache.json').read_text())
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
assert subprocess.check_output(aws + ['sts', 'get-caller-identity', '--query', 'Account', '--output', 'text'], text=True).strip() == '839993661825'
base = 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
prefix = base + 'protein-nr-query-calls/'
subprocess.run(aws + ['s3', 'cp', str(cache) + '/', prefix, '--recursive', '--no-follow-symlinks', '--only-show-errors'], check=True)
subprocess.run(aws + ['s3', 'cp', str(d / 'query-call-cache.json'), prefix + 'query-call-cache.json', '--only-show-errors'], check=True)
subprocess.run(aws + ['s3', 'cp', str(r / 'protein/nr-query-call-full.time.txt'), prefix + 'time.txt', '--only-show-errors'], check=True)
checks = [('raw_NR', base + 'protein-nr-combined/nr.tsv', m['raw'])]
checks += [(kind, prefix + Path(x['path']).name, x) for kind, x in m['files'].items()]
result = {'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'account': '839993661825', 'method': 'Stream each S3 object through SHA256 and compare byte count and digest with the completed local source manifest.', 'files': []}
for kind, uri, expected in checks:
    start = time.monotonic()
    p = subprocess.Popen(aws + ['s3', 'cp', uri, '-', '--only-show-errors'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    h, size = hashlib.sha256(), 0
    for block in iter(lambda: p.stdout.read(16 << 20), b''):
        h.update(block); size += len(block)
    error = p.stderr.read().decode()
    assert p.wait() == 0, error
    assert size == expected['bytes'] and h.hexdigest() == expected['sha256'], f'Download mismatch: {kind}'
    item = dict(kind=kind, uri=uri, bytes=size, sha256=h.hexdigest(), seconds=time.monotonic() - start)
    result['files'].append(item)
    print(json.dumps(item), flush=True)
    (cache / 'download-verification-progress.json').write_text(json.dumps(result, indent=2) + '\n')
result['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
result['all_download_hashes_verified'] = True
receipt = cache / 'backup-download-verification.json'
receipt.write_text(json.dumps(result, indent=2) + '\n')
subprocess.run(aws + ['s3', 'cp', str(receipt), prefix + receipt.name, '--only-show-errors'], check=True)
print('NR_RAW_AND_COMPACT_BACKUPS_DOWNLOAD_HASH_VERIFIED', flush=True)

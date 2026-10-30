"""Upload and download-hash-verify a completed pilot snapshot to private S3."""
import argparse
import datetime
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--directory', type=Path, default=Path('/home/ec2-user/experiment/cohort-storage/translated-core-nt/pilot24'))
p.add_argument('--label', required=True)
p.add_argument('--commit', required=True)
a = p.parse_args()
assert re.fullmatch('[a-z0-9-]+', a.label) and re.fullmatch('[a-f0-9]{40}', a.commit)
assert a.directory.name == 'pilot24'
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
assert subprocess.check_output(aws + ['sts', 'get-caller-identity', '--query', 'Account', '--output', 'text'], text=True).strip() == '839993661825'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
for method in ['megablast', 'blastn', 'tblastx']:
    d = a.directory / method
    c = json.loads((d / 'complete.json').read_text())
    assert json.loads((d / 'resources.json').read_text())['status'] == 'complete'
    assert sha(d / 'alignments.tsv') == c['output_sha256']
prefix = 's3://pipeline-lab-839993661825-us-east-1/runs/translated-core-nt-20261022/pilot-snapshots/' + a.label + '/'
bucket, key_prefix = prefix[5:].split('/', 1)
assert subprocess.check_output(aws + ['s3api', 'list-objects-v2', '--bucket', bucket, '--prefix', key_prefix, '--max-keys', '1', '--query', 'KeyCount', '--output', 'text'], text=True).strip() == '0', 'Snapshot label already used; choose a fresh label'
receipt_path = a.directory.parent / ('backup-' + a.label + '.json')
files = []
with tempfile.TemporaryDirectory() as tmp:
    for path in sorted(a.directory.rglob('*')):
        assert not path.is_symlink(), path
        if not path.is_file():
            continue
        relative = str(path.relative_to(a.directory))
        digest = sha(path)
        subprocess.run(aws + ['s3', 'cp', str(path), prefix + relative, '--only-show-errors'], check=True)
        downloaded = Path(tmp) / 'check'
        subprocess.run(aws + ['s3', 'cp', prefix + relative, str(downloaded), '--only-show-errors'], check=True)
        assert sha(downloaded) == digest and sha(path) == digest
        files.append(dict(path=relative, bytes=path.stat().st_size, sha256=digest))
receipt = dict(verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_commit=a.commit, backup_script_sha256=sha(Path(__file__)), prefix=prefix, analysis_present=(a.directory / 'analysis/summary.json').exists(), files=files, verification='Every selected file uploaded, downloaded and SHA256-compared; source unchanged during verification. No additional search launched.')
receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
subprocess.run(aws + ['s3', 'cp', str(receipt_path), prefix + 'backup-verification.json', '--only-show-errors'], check=True)
print(json.dumps({k: v for k, v in receipt.items() if k != 'files'}), flush=True)

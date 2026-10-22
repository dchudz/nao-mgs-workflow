"""Preserve completed primary-worker evidence without archiving transient task work directories."""
import datetime, hashlib, json, os, subprocess, tarfile
from pathlib import Path

r = Path('/home/ec2-user/experiment')
closure = r / 'primary-closure'
closure.mkdir(exist_ok=True)
base = 'runs/amino-acid-20261021/'
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
assert subprocess.check_output(aws + ['sts', 'get-caller-identity', '--query', 'Account', '--output', 'text'], text=True).strip() == '839993661825'
scopes = [
    ('cohort/heldout', 'cohort-screens'),
    ('cohort/development', 'primary-development-evidence'),
    ('full/comparison', 'full-depth-comparison'),
    ('full/protein', 'full-protein'),
    *[('full/' + arm, 'full-depth-evidence/' + arm) for arm in ['control', 'aa12', 'nt20', 'nt22']],
    ('pilot', 'pilot-evidence'),
    ('null-controls', 'null-controls-evidence'),
    ('cohort-storage/protein200', 'protein200'),
    ('cohort-storage/nr-native-controls', 'nr-native-controls'),
    ('cohort-storage/seed-diagnostic200', 'seed-diagnostic200'),
    ('cohort-storage/validation-scaling', 'validation-scaling'),
    ('cohort-storage/repeat-variant', 'repeat-variant'),
    ('cohort-storage/aa-loss-diagnostic', 'aa-loss-diagnostic'),
    ('cohort-storage/aa-over-nt20', 'aa-over-nt20'),
    ('cohort-storage/candidate-outlier', 'candidate-outlier'),
    ('cohort-storage/screen-burden-inputs', 'screen-burden-inputs'),
]
records = []
for relative, prefix in scopes:
    root = r / relative
    assert root.is_dir(), root
    exclusions = []
    for directory, names, files in os.walk(root):
        for name in list(names):
            if name == 'work' or name.endswith('-work'):
                exclusions.append(str((Path(directory) / name).relative_to(root)))
                names.remove(name)
    sync = aws + ['s3', 'sync', str(root) + '/', 's3://pipeline-lab-839993661825-us-east-1/' + base + prefix + '/', '--size-only', '--only-show-errors']
    verify = ['python3', str(r / 'verify-private-prefix.py'), str(root), base + prefix]
    for name in exclusions:
        sync.extend(['--exclude', name + '/*'])
        verify.extend(['--exclude-dir', name])
    result = closure / ('static-' + relative.replace('/', '-') + '.json')
    verify.extend(['--output', str(result)])
    subprocess.run(sync, check=True)
    subprocess.run(verify, check=True, stdout=subprocess.DEVNULL)
    record = json.loads(result.read_text())
    records.append({k: record[k] for k in ['local_root', 'prefix', 'expected_files', 'expected_bytes', 'verified']})
    print(json.dumps(records[-1]), flush=True)

# Root scripts, parameter records and logs are compact; public reference caches and
# downloaded raw runs remain reconstructible from their preserved source manifests.
files = [p for p in r.iterdir() if p.is_file() and p.suffix in ['.py', '.json', '.jsonl', '.log', '.tsv', '.nf', '.sh', '.rs', '.txt']]
manifest = {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
(closure / 'primary-root-manifest.json').write_text(json.dumps(manifest, indent=2))
archive = closure / 'primary-root-evidence.tar.gz'
with tarfile.open(archive, 'w:gz') as tar:
    for p in files:
        tar.add(p, arcname=p.name, recursive=False)
record = {'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'scopes': records, 'archive_bytes': archive.stat().st_size, 'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'scope': 'Completed static scientific evidence. Mutable protein final analyses and native-validation results require their separate completion and backup checks before worker retirement.'}
(closure / 'static-backups-complete.json').write_text(json.dumps(record, indent=2))
subprocess.run(aws + ['s3', 'cp', str(closure) + '/', 's3://pipeline-lab-839993661825-us-east-1/' + base + 'primary-closure/', '--recursive', '--only-show-errors'], check=True)
print('PRIMARY_STATIC_BACKUPS_VERIFIED', flush=True)

from pathlib import Path
import datetime, hashlib, json, os, subprocess, tarfile

root = Path('/home/ec2-user/experiment')
out = root / 'preservation-20261024'
out.mkdir(exist_ok=True)
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
identity = json.loads(subprocess.check_output(aws + ['sts', 'get-caller-identity']))
assert identity['Account'] == '839993661825'
prefix = 's3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/preservation/20261024/'
sample = root / 'e2e/SRR37076045'
work = sample / 'protein-run/work/62/f6a71d0c492f2c3e0e47e2f547717e'
assert (work / '.exitcode').read_text().strip() == '0'
assert not subprocess.check_output(['docker', 'ps', '-q']).strip()

def digest_file(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()

records = []
def preserve(path, name):
    before = path.stat()
    sha = digest_file(path)
    uri = prefix + name
    print(json.dumps({'phase': 'upload', 'name': name, 'bytes': before.st_size}), flush=True)
    subprocess.run(aws + ['s3', 'cp', str(path), uri, '--only-show-errors'], check=True)
    print(json.dumps({'phase': 'independent_readback', 'name': name}), flush=True)
    proc = subprocess.Popen(aws + ['s3', 'cp', uri, '-', '--only-show-errors'], stdout=subprocess.PIPE)
    remote = hashlib.sha256()
    size = 0
    for block in iter(lambda: proc.stdout.read(16 * 1024**2), b''):
        size += len(block)
        remote.update(block)
    proc.stdout.close()
    assert proc.wait() == 0
    assert size == before.st_size and remote.hexdigest() == sha
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    records.append({'uri': uri, 'bytes': size, 'sha256': sha, 'readback_verified': True})
    (out / 'backup-verification.json').write_text(json.dumps(records, indent=2) + '\n')
    print(json.dumps({'phase': 'verified', **records[-1]}), flush=True)

archive = out / 'completed-artifacts-and-diagnostics.tar.gz'
with tarfile.open(archive, 'w:gz', compresslevel=1) as tar:
    for directory, dirs, names in os.walk(root / 'e2e'):
        dirs[:] = [name for name in dirs if name not in {'work', '.nextflow', 'checkpoint-staging'}]
        for name in names:
            path = Path(directory) / name
            if path.is_file():
                tar.add(path, arcname=str(path.relative_to(root)), recursive=False)
    for directory in (root / 'e2e').glob('*/*/work'):
        for path in directory.glob('*/*/*'):
            if path.is_file() and not path.is_symlink() and (path.name.startswith('.command.') or path.name == '.exitcode' or path.suffix == '.json' or path.name.endswith('.stats.txt')):
                tar.add(path, arcname=str(path.relative_to(root)), recursive=False)
    for path in [root / '.nextflow', root / 'experiments/amino_acids/end-to-end-full-depth', root / 'setup/progress-backup-latest.json']:
        if path.exists():
            tar.add(path, arcname=str(path.relative_to(root)))
    for path in (root / 'data').glob('*/verified.json'):
        tar.add(path, arcname=str(path.relative_to(root)))
preserve(archive, archive.name)
preserve(work / 'protein-hits.tsv', 'completed-sorted-search-output.tsv')
receipt = {'verified_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'records': records,
           'all_four_outputs_completed': False, 'raw_inputs_and_remaining_scratch': 'Retained on attached EBS volume; not duplicated in this archive.',
           'deletions': [], 'new_searches_launched': False}
(out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
subprocess.run(aws + ['s3', 'cp', str(out / 'receipt.json'), prefix + 'receipt.json', '--only-show-errors'], check=True)
print(json.dumps({'phase': 'complete', **receipt}), flush=True)

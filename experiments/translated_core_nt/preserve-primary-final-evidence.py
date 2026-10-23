"""Archive completed primary-worker evidence; never stop or delete resources.

Run only after scientific analyses finish. Each archive records every source hash
and is independently downloaded and hashed before its receipt is accepted.
"""
import datetime
import hashlib
import json
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path

r = Path('/home/ec2-user/experiment')
d = r / 'cohort-storage/translated-core-nt'
out = r / 'primary-closure-final'
aws = ['aws', '--profile', 'pipeline-codex', '--region', 'us-east-1']
prefix = 's3://pipeline-lab-839993661825-us-east-1/runs/translated-core-nt-20261022/primary-closure-final/'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def verify_upload(path):
    digest = sha(path)
    subprocess.run(aws + ['s3', 'cp', str(path), prefix + path.name, '--only-show-errors'], check=True)
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        copy = Path(tmp) / path.name
        subprocess.run(aws + ['s3', 'cp', prefix + path.name, str(copy), '--only-show-errors'], check=True)
        assert copy.stat().st_size == path.stat().st_size and sha(copy) == digest
    return dict(uri=prefix + path.name, bytes=path.stat().st_size,
                sha256=digest, independent_download_hash_verified=True)


assert subprocess.check_output(aws + ['sts', 'get-caller-identity', '--query', 'Account', '--output', 'text'], text=True).strip() == '839993661825'
for process in ['blastn', 'tblastx', 'tblastn', 'diamond']:
    assert not subprocess.run(['pgrep', '-x', process], capture_output=True).stdout.strip(), process
for method in ['megablast', 'blastn', 'tblastx-word4']:
    resource = json.loads((d / 'protein232' / method / 'resources.json').read_text())
    assert resource['status'] in ['complete', 'failed', 'time_limit'], resource
assert (r / 'cohort-storage/protein-finalization/complete.json').exists()
core = json.loads((d / 'reference-archive/complete.json').read_text())
assert core['verified'] and len(core['files']) == 836
out.mkdir(exist_ok=True)
assert not (out / 'complete.json').exists(), 'A final closure already exists'

metadata_dirs = ['.nextflow', 'assets', 'bin', 'configs', 'experiments', 'lib',
                 'modules', 'schemas', 'subworkflows', 'workflows', 'rust-tools',
                 'tools', 'ref', 'protein-repeat-gate', 'protein-screen-boundaries',
                 'range-probe', 'reads', 'reads-pilot', 'primary-closure', 'backup-verification']
scopes = {name.replace('/', '-'): [r / name] for name in
          ['cohort-storage', 'cohort', 'full', 'pilot', 'protein', 'reference/results', 'null-controls']}
scopes['metadata-and-inputs'] = [r / name for name in metadata_dirs] + [p for p in r.iterdir() if p.is_file()]
skip = {r / 'reference/results/blast_db': 'Separately checksum-verified full core_nt Deep Archive snapshot',
        r / 'protein/prot.accession2taxid.FULL.gz': 'Unused downloaded accession mapping cache; original download log/source retained'}
records = []
for name, roots in scopes.items():
    receipt = out / (name + '.json')
    target = out / (name + '.tar.gz')
    if receipt.exists():
        previous = json.loads(receipt.read_text())
        assert previous['archive']['independent_download_hash_verified']
        assert target.stat().st_size == previous['archive']['bytes'] and sha(target) == previous['archive']['sha256']
        records.append(previous)
        continue
    assert not target.exists(), ('Preserve interrupted archive before retrying', target)
    paths, exclusions = set(), []
    for root in roots:
        assert root.exists(), root
        if root.is_file() or root.is_symlink():
            paths.add(root)
            continue
        for directory, dirs, files in os.walk(root, followlinks=False):
            for child in list(dirs):
                p = Path(directory) / child
                if p in skip or child == 'work' or child.endswith('-work') or child == '__pycache__':
                    dirs.remove(child)
                    exclusions.append(dict(path=str(p.relative_to(r)), reason=skip.get(p, 'Transient workflow work or Python cache')))
                elif p.is_symlink():
                    dirs.remove(child)
                    paths.add(p)
            for child in files:
                p = Path(directory) / child
                if p in skip:
                    exclusions.append(dict(path=str(p.relative_to(r)), reason=skip[p]))
                else:
                    paths.add(p)
    manifest = []
    class Reader:
        def __init__(self, f):
            self.f, self.digest, self.bytes = f, hashlib.sha256(), 0
        def read(self, n):
            b = self.f.read(n)
            self.digest.update(b)
            self.bytes += len(b)
            return b
    with tarfile.open(target, 'w:gz', compresslevel=1) as tar:
        for path in sorted(paths):
            relative = str(path.relative_to(r))
            if path.is_symlink():
                tar.add(path, arcname=relative, recursive=False)
                manifest.append(dict(path=relative, symlink=os.readlink(path),
                                     resolved_target=str(path.resolve()), target_exists=path.exists()))
                continue
            before = path.stat()
            info = tar.gettarinfo(str(path), arcname=relative)
            with path.open('rb') as f:
                reader = Reader(f)
                tar.addfile(info, reader)
            assert reader.bytes == before.st_size
            # Only continuously appended root logs are intentionally snapshots.
            root_log = path.parent == r and path.suffix in ('.log', '.jsonl')
            after = path.stat()
            if not root_log:
                assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), path
            manifest.append(dict(path=relative, bytes=reader.bytes,
                                 sha256=reader.digest.hexdigest(), root_log_snapshot=root_log))
    record = dict(scope=name, completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  archive=verify_upload(target), files=manifest, exclusions=exclusions)
    receipt.write_text(json.dumps(record, indent=2) + '\n')
    verify_upload(receipt)
    records.append(record)
    print(json.dumps({'completed_scope': name, 'archive': record['archive'], 'source_files': len(manifest)}), flush=True)

final = dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             account='839993661825', source_sha256=sha(Path(__file__)),
             archives=[dict(scope=x['scope'], **x['archive']) for x in records],
             scope_receipts=[x['scope'] + '.json' for x in records],
             frozen_core_reference_prefix=core['prefix'],
             all_scientific_archives_download_hash_verified=True,
             limits=['Transient work directories and public SRA download caches are excluded.',
                     'Reference blast_db is preserved separately with recorded SHA and end-to-end CRC verification; no archived-reference restore was performed.',
                     'Existing S3 scientific backups are unchanged. This closure adds consolidated source-hashed archives.',
                     'Symlinks are retained as links with their original target recorded. Inspect their receipts before resource retirement.',
                     'No resources are stopped, deleted or made public.'])
(out / 'complete.json').write_text(json.dumps(final, indent=2) + '\n')
verify_upload(out / 'complete.json')
print(json.dumps(final), flush=True)

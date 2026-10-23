"""Write source-hashed tar entries, retaining symlink and hardlink semantics."""
import hashlib
import os
import subprocess
import tarfile
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


class Reader:
    def __init__(self, stream):
        self.stream, self.digest, self.bytes = stream, hashlib.sha256(), 0

    def read(self, n):
        data = self.stream.read(n)
        self.digest.update(data)
        self.bytes += len(data)
        return data


def archive_paths(target, paths, root):
    manifest, added = [], {}
    with target.open('xb') as sink:
        process = subprocess.Popen(['/usr/bin/pigz', '-n', '-1', '-p', '8'],
                                   stdin=subprocess.PIPE, stdout=sink)
        try:
            with tarfile.open(fileobj=process.stdin, mode='w|') as tar:
                for path in sorted(paths):
                    relative = str(path.relative_to(root))
                    if path.is_symlink():
                        tar.add(path, arcname=relative, recursive=False)
                        manifest.append(dict(path=relative, symlink=os.readlink(path),
                                             resolved_target=str(path.resolve()), target_exists=path.exists()))
                        continue
                    before = path.stat()
                    info = tar.gettarinfo(str(path), arcname=relative)
                    if info.islnk():
                        # A tar hardlink stores a target name and no second payload.
                        sha = digest(path)
                        assert added[info.linkname]['sha256'] == sha
                        assert added[info.linkname]['bytes'] == before.st_size
                        tar.addfile(info)
                        row = dict(path=relative, bytes=before.st_size, sha256=sha, hardlink=info.linkname)
                    else:
                        assert info.isfile(), path
                        with path.open('rb') as stream:
                            reader = Reader(stream)
                            tar.addfile(info, reader)
                        assert reader.bytes == before.st_size
                        row = dict(path=relative, bytes=reader.bytes, sha256=reader.digest.hexdigest())
                    root_log = path.parent == root and path.suffix in ('.log', '.jsonl')
                    after = path.stat()
                    if not root_log:
                        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), path
                    row['root_log_snapshot'] = root_log
                    manifest.append(row)
                    added[relative] = row
        finally:
            try:
                process.stdin.close()
            finally:
                status = process.wait()
        assert status == 0, ('Parallel gzip failed', status)
    return manifest

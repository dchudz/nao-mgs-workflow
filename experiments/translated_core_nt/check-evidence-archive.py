"""A restored archive retains bytes and link identity for every fixture entry."""
import hashlib,json,os,subprocess,tarfile,tempfile
from pathlib import Path
from evidence_archive import archive_paths

with tempfile.TemporaryDirectory() as tmp:
    base=Path(tmp);source=base/'source';source.mkdir();(source/'nested').mkdir()
    data=bytes(range(256))*4096
    (source/'a.bin').write_bytes(data);os.link(source/'a.bin',source/'nested/b.bin')
    (source/'empty').touch();(source/'reference').symlink_to('a.bin')
    paths={source/'a.bin',source/'nested/b.bin',source/'empty',source/'reference'}
    original={str(p.relative_to(source)):p.read_bytes() for p in paths if not p.is_symlink()}
    archive=base/'evidence.tar.gz';manifest=archive_paths(archive,paths,source)
    with tarfile.open(archive,'r:gz') as tar:
        members={x.name:x for x in tar.getmembers()}
        assert members['nested/b.bin'].islnk() and members['nested/b.bin'].size==0
        assert members['reference'].issym()
        restored=base/'restored';restored.mkdir()
        # Only this locally generated fixture is extracted.
        tar.extractall(restored)
    for name,expected in original.items():assert (restored/name).read_bytes()==expected
    assert (restored/'a.bin').stat().st_ino==(restored/'nested/b.bin').stat().st_ino
    assert (restored/'reference').is_symlink() and os.readlink(restored/'reference')=='a.bin'
    for row in manifest:
        if 'sha256' in row:assert row['sha256']==hashlib.sha256((restored/row['path']).read_bytes()).hexdigest()
    assert all(p.read_bytes()==original[str(p.relative_to(source))] for p in paths if not p.is_symlink())
    result={'regular_file_bytes_preserved':True,'empty_file_preserved':True,'hardlink_identity_preserved':True,'symlink_identity_preserved':True,'every_manifest_hash_matches_restored_bytes':True,'source_contents_unchanged':True,'hardlink_tar_payload_bytes':members['nested/b.bin'].size,'hardlink_source_bytes':len(data),'compression_command':['/usr/bin/pigz','-n','-1','-p','8'],'compressor_version':subprocess.check_output(['/usr/bin/pigz','--version'],text=True).strip(),'helper_sha256':hashlib.sha256(Path(__file__).with_name('evidence_archive.py').read_bytes()).hexdigest()}
    print(json.dumps(result))

"""Locate a numeric DIAMOND target in a flat native BLAST protein database.

Reads PIN headers only. It does not search, decode sequences, or change the database.
The field layout and cumulative OID convention follow pinned DIAMOND v2.1.17
src/data/blastdb/{pin,pal}.cpp. Unsupported filtered/nested aliases fail closed.
"""
import argparse
import bisect
import hashlib
import json
import shlex
import struct
from pathlib import Path


def volume_header(path):
    with path.open('rb') as f:
        def exact(n):
            b = f.read(n)
            assert len(b) == n, (path, 'truncated header')
            return b
        def u32():
            return struct.unpack('>I', exact(4))[0]
        def string():
            n = u32()
            assert n < 1024**2, (path, 'unbounded header string')
            return exact(n).decode('utf-8').rstrip('\0')
        version, molecule = u32(), u32()
        assert version in (4, 5) and molecule == 1, (path, version, molecule)
        number = u32() if version == 5 else None
        title = string()
        lmdb = string() if version == 5 else None
        date = string()
        count = u32()
        letters = struct.unpack('<Q', exact(8))[0]
        maximum_length = u32()
        header_bytes = f.tell()
        assert count > 0 and path.stat().st_size >= header_bytes + 8 * (count + 1)
        f.seek(0)
        digest = hashlib.sha256(f.read(header_bytes)).hexdigest()
    return dict(version=version, volume_number=number, sequences=count,
                letters=letters, maximum_length=maximum_length, title=title,
                date=date, lmdb=lmdb, header_bytes=header_bytes,
                header_sha256=digest, pin_bytes=path.stat().st_size)


def locate(db, ordinals):
    alias = Path(str(db) + '.pal')
    if alias.exists():
        fields = {}
        for line in alias.read_text().splitlines():
            values = shlex.split(line, comments=True)
            if not values:
                continue
            assert values[0] in ('TITLE', 'DBLIST', 'NSEQ', 'LENGTH'), values[0]
            assert values[0] not in fields, values[0]
            fields[values[0]] = values[1:]
        names = fields['DBLIST']
        assert names and len(names) == len(set(names))
        volumes = [alias.parent / name for name in names]
        assert not any(Path(str(v) + '.pal').exists() for v in volumes), 'Nested alias'
    else:
        fields, volumes = {}, [db]
    records, total, letters = [], 0, 0
    for v in volumes:
        record = volume_header(Path(str(v) + '.pin'))
        record.update(database=str(v), first_global_ordinal=total)
        total += record['sequences']
        letters += record['letters']
        record['end_global_ordinal_exclusive'] = total
        records.append(record)
    for key, actual in [('NSEQ', total), ('LENGTH', letters)]:
        if key in fields:
            assert len(fields[key]) == 1 and int(fields[key][0]) == actual
    starts = [x['first_global_ordinal'] for x in records]
    targets = []
    for oid in ordinals:
        assert 0 <= oid < total, oid
        v = records[bisect.bisect_right(starts, oid) - 1]
        local = oid - v['first_global_ordinal']
        targets.append(dict(global_ordinal=oid, database=v['database'],
                            local_ordinal=local, diamond_getseq_number=local + 1))
    return dict(database=str(db), sequences=total, letters=letters,
                volumes=records, targets=targets,
                rule='snum is zero-based; getseq --seq is one-based within the opened database. A single-volume getseq avoids scanning all NR volumes. This is retrieval after the comprehensive search, not restricted-reference validation.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('database', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--ordinals', type=int, nargs='*', default=[])
    a = p.parse_args()
    result = locate(a.database, a.ordinals)
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: result[k] for k in ('sequences', 'letters', 'targets')}))

"""Read only verified, equivalent NR evidence caches; retain raw-input fallback."""
import csv
import hashlib
import json
from pathlib import Path

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''): h.update(block)
    return h.hexdigest()

def manifest(directory):
    d = Path(directory)
    p = d / 'query-call-cache.json'
    if not p.exists(): return None
    m = json.loads(p.read_text())
    assert m['format'] == 'NR near-best query calls v1'
    raw = (d / 'nr.tsv').stat()
    assert raw.st_size == m['raw']['bytes'] and raw.st_mtime_ns == m['raw']['mtime_ns'], 'NR source changed after compaction'
    assert m['parity']['all_query_taxa_scores_and_classifications_identical']
    return m

def verified_file(m, key):
    x = m['files'][key]
    p = Path(x['path'])
    assert p.stat().st_size == x['bytes'] and digest(p) == x['sha256'], 'Evidence cache hash mismatch'
    return p

def alignment_input(directory, kind):
    assert kind in ['native_controls', 'seeds']
    m = manifest(directory)
    return verified_file(m, kind) if m else Path(directory) / 'nr.tsv'

def near_best_calls(directory, original):
    m = manifest(directory)
    if not m:
        with (Path(directory) / 'nr.tsv').open() as f:
            return original.near_best_query_taxa(csv.reader(f, delimiter='\t'))
    calls = {}
    with verified_file(m, 'query_calls').open() as f:
        for line in f:
            q, score, taxa = line.rstrip('\n').split('\t')
            assert q not in calls
            values = {int(t) for t in taxa.split(';')}
            assert values
            calls[q] = (values, float(score))
    assert len(calls) == m['counts']['near_best_query_rows']
    return calls

"""Recover original candidate-read qualities with complete membership/sequence checks."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, default=Path('/home/ec2-user/experiment'))
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
r, out = a.root, a.output
out.mkdir(parents=True, exist_ok=False)
sys.path.insert(0, str(r / 'experiments/amino_acids'))
from union_pairs import pairs
def table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        yield from csv.DictReader(f, delimiter='\t')
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()
d = r / 'cohort-storage/protein200'
viral = {x['query_pair'] for x in table(d / 'viral-pairs.tsv')}
meta = {x['query_pair']: x for x in table(d / 'all.tsv.gz') if x['query_pair'] in viral}
assert meta.keys() == viral
members = collections.defaultdict(dict)
for x in table(d / 'all.members.tsv.gz'):
    q = x['query_pair']
    if q in viral:
        sample = q.split('_', 1)[0]
        assert x['seq_id'] not in members[sample]
        members[sample][x['seq_id']] = q
samples = [x['run_accession'] for x in json.loads((r / 'heldout-ranges.json').read_text())][:200]
assert len(samples) == len(set(samples)) == 200 and members.keys() <= set(samples)
batches = {s: p.parent for p in (r / 'cohort-storage/batches').glob('*/complete.json') for s in json.loads(p.read_text())['samples']}
metrics, total = [], 0
with gzip.open(out / 'qualities.tsv.gz', 'wt') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['query_pair', 'query_qual', 'query_qual_rev', 'source_seq_id', 'observed_original_pairs', 'sequence_pair_sha256'])
    for sample in samples:
        start = time.monotonic()
        source = batches[sample] / 'output/fastp' / f'{sample}_aa12_fastp.fastq.gz'
        wanted, seen, best, observations = members[sample], set(), {}, collections.Counter()
        for key, _ in pairs(source):
            name, s1, q1, s2, q2 = [x.decode() for x in key]
            seq_id = name.lstrip('@')
            if seq_id not in wanted:
                continue
            assert seq_id not in seen
            seen.add(seq_id)
            q = wanted[seq_id]
            assert (s1, s2) == (meta[q]['query_seq'], meta[q]['query_seq_rev']), (sample, seq_id, 'Sequence mismatch')
            candidate = sum(map(ord, q1 + q2)), q1, q2, name
            if q not in best or candidate > best[q]:
                best[q] = candidate
            observations[q] += 1
        assert seen == wanted.keys(), (sample, 'Original FASTP membership mismatch')
        for q in sorted(best):
            assert observations[q] == int(meta[q]['occurrences'])
            _, q1, q2, name = best[q]
            digest = hashlib.sha256((meta[q]['query_seq'] + '\n' + meta[q]['query_seq_rev']).encode()).hexdigest()
            writer.writerow([q, q1, q2, name, observations[q], digest])
        total += len(best)
        metric = dict(sample=sample, candidate_pairs=len(best), original_pairs=len(seen), source=str(source), source_sha256=sha(source), seconds=time.monotonic() - start)
        metrics.append(metric)
        print(json.dumps(metric), flush=True)
assert total == len(viral)
result = dict(candidate_pairs=total, original_pairs=sum(x['original_pairs'] for x in metrics), samples=len(samples), source_sha256=sha(Path(__file__)), metadata_sha256=sha(d / 'all.tsv.gz'), membership_sha256=sha(d / 'all.members.tsv.gz'), viral_calls_sha256=sha(d / 'viral-pairs.tsv'), qualities_sha256=sha(out / 'qualities.tsv.gz'), sample_metrics=metrics, policy='Every candidate member must occur once with exactly matching sequences in the original AA12 FASTP source. Select maximum summed original qualities per exact pair, with lexical quality and read-ID tie breaking. No rerun, synthetic qualities, or taxonomic verdict changes.')
(out / 'complete.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'sample_metrics'}))

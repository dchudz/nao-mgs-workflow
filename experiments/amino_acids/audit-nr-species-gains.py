"""Inspect fixed-hash examples of species gains using existing viral alignments only."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


def table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        yield from csv.DictReader(f, delimiter='\t')


def ranked(rows, key, salt):
    return sorted(rows, key=lambda x: hashlib.sha256((salt + key(x)).encode()).digest())


def coordinates(subject, start, end, genome_length):
    genome, frame, offset = subject.rsplit('~', 2)
    frame, offset = int(frame), int(offset)
    assert 0 <= frame < 6 and 1 <= start <= end
    lo, hi = offset + 3 * (start - 1), offset + 3 * end
    assert 0 <= lo < hi <= genome_length
    if frame >= 3:
        lo, hi = genome_length - hi, genome_length - lo
    return genome, '-' if frame >= 3 else '+', lo, hi


def main(a):
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    groups = collections.defaultdict(list)
    for x in table(a.calls):
        if x['baseline_arm'] != 'nt20' or x['tier'] != 'NR_and_viral_reference_same_species' or x['gained_protein_call'] != 'True':
            continue
        group = 'zero_baseline_support' if x['baseline_exact_pairs'] == '0' else 'crosses_three_pair_threshold'
        groups[group].append(x)
    selected = {}
    for group, rows in sorted(groups.items()):
        for x in ranked(rows, lambda x: x['sample'] + ':' + x['species'], 'species-gain-call:')[:5]:
            selected[x['sample'], x['species']] = dict(x, audit_group=group, stratum_calls=len(rows))
    candidates = collections.defaultdict(list)
    for x in table(a.resolution):
        key = x['sample'], x['nr_species']
        if key in selected and x['species_relation'] == 'same_species':
            candidates[key].append(x)
    queries = {}
    for key, call in selected.items():
        rows = candidates[key]
        assert len(rows) == int(call['protein_exact_pairs']), (key, len(rows), call)
        for x in ranked(rows, lambda x: x['query_pair'], 'species-gain-pair:')[:12]:
            queries[x['query_pair']] = dict(x, audit_group=call['audit_group'])
    # Keep five deterministic highest-scoring rows per mate. Count all near-best
    # rows in a second pass over the selected-query extraction, not the 27 GB input.
    best = collections.defaultdict(float)
    extracted = out / 'selected-viral-alignments.tsv'
    with open(a.hits) as source, extracted.open('w') as dest:
        for line in (source if queries else ()):
            q = line.split('\t', 1)[0]
            if q.rsplit('/', 1)[0] in queries:
                dest.write(line)
                row = line.rstrip('\n').split('\t')
                assert len(row) == 12
                if int(row[4]) >= 25:
                    best[q] = max(best[q], float(row[2]))
    kept = collections.defaultdict(list)
    near_counts = collections.Counter()
    with extracted.open() as f:
        for row in csv.reader(f, delimiter='\t'):
            q = row[0]
            if int(row[4]) < 25 or float(row[2]) < .9 * best[q]:
                continue
            near_counts[q] += 1
            kept[q].append(row)
            kept[q].sort(key=lambda x: (-float(x[2]), x[1], int(x[9]), int(x[10]), int(x[7]), int(x[8])))
            del kept[q][5:]
    wanted = {row[1].rsplit('~', 2)[0] for rows in kept.values() for row in rows}
    lengths = {}
    name, n = None, 0
    with gzip.open(a.genomes, 'rt') as f:
        for line in (f if wanted else ()):
            if line.startswith('>'):
                if name in wanted:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
        if name in wanted:
            lengths[name] = n
    assert wanted == lengths.keys()
    hits = []
    for q, rows in sorted(kept.items()):
        for rank, row in enumerate(rows, 1):
            genome = row[1].rsplit('~', 2)[0]
            genome, strand, start, end = coordinates(row[1], int(row[9]), int(row[10]), lengths[genome])
            hits.append(dict(query=q, retained_rank=rank, eligible_near_best_rows=near_counts[q], subject=row[1], genome=genome, reference_strand=strand, reference_start_0based=start, reference_end_exclusive=end, bitscore=row[2], evalue=row[3], alignment_aa=row[4], identity_pct=row[5], query_length_nt=row[6], query_start=row[7], query_end=row[8], query_frame=row[11]))
    ids = collections.defaultdict(list)
    for x in (table(a.members) if queries else ()):
        if x['query_pair'] in queries:
            ids[x['query_pair']].append(x['seq_id'])
            ids[x['query_pair']] = sorted(ids[x['query_pair']])[:3]
    for q, x in queries.items():
        assert ids[q] and any(q + '/' + str(m) in kept for m in [1, 2])
        x['example_original_read_ids'] = ';'.join(ids[q])
    for name, rows in [('calls.tsv', list(selected.values())), ('pairs.tsv', list(queries.values())), ('alignments.tsv', hits)]:
        with (out / name).open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ['none'], delimiter='\t')
            w.writeheader()
            w.writerows(rows)
    summary = dict(selected_calls=len(selected), selected_pairs=len(queries), retained_alignment_rows=len(hits), source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), call_table_sha256=hashlib.sha256(Path(a.calls).read_bytes()).hexdigest(), resolution_sha256=hashlib.sha256(Path(a.resolution).read_bytes()).hexdigest(), selection='Up to five fixed-hash NT20-relative gained protein calls per baseline-support stratum, restricted to concordant NR and viral-reference species; up to twelve fixed-hash pairs per call.', limits=['Illustrative outcome-selected cases, not a representative accuracy estimate.', 'Three exact pairs is not the production alignment/similarity exemplar threshold.', 'Reference intervals come from existing viral-reference alignments, not NR alignments.', 'Up to five highest-scoring eligible rows retained per mate; ties use subject and coordinates deterministically. The selected-query raw extraction retains every row.', 'Overlapping intervals can identify concentrated support but cannot prove molecule independence or infection.'])
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['calls', 'resolution', 'hits', 'genomes', 'members', 'output']:
        p.add_argument('--' + name, required=True)
    main(p.parse_args())

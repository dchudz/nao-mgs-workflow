"""Check raw output multiplicity and every frozen query for a batching replay."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

FIELDS = 'qseqid sseqid staxids bitscore evalue length pident qcovs qcovhsp qlen slen qstart qend sstart send qframe sframe qseq sseq'.split()


def compare(left, right, queries):
    tables = []
    for path in (left, right):
        rows = collections.Counter()
        with path.open() as f:
            for row in csv.reader(f, delimiter='\t'):
                assert len(row) == len(FIELDS) and row[0] in queries, (path, row)
                rows[tuple(row)] += 1
        tables.append(rows)
    a, b = tables
    only = (a-b, b-a)
    counts = [collections.Counter(), collections.Counter()]
    differences = [collections.Counter(), collections.Counter()]
    for side in range(2):
        for row, n in tables[side].items():
            counts[side][row[0]] += n
        for row, n in only[side].items():
            differences[side][row[0]] += n
    per_query = [dict(query=q, left_rows=counts[0][q], right_rows=counts[1][q],
                      only_left_rows=differences[0][q], only_right_rows=differences[1][q],
                      multiset_equal=not (differences[0][q] or differences[1][q]))
                 for q in queries]
    return dict(alignment_multiset_equal=a == b, query_count=len(queries),
                left_rows=sum(a.values()), right_rows=sum(b.values()),
                rows_only_left=sum(only[0].values()), rows_only_right=sum(only[1].values()),
                differing_queries=sum(not x['multiset_equal'] for x in per_query),
                no_output_queries_left=sum(not x['left_rows'] for x in per_query),
                no_output_queries_right=sum(not x['right_rows'] for x in per_query)), per_query, only


def main(directory, left, right, output):
    raw = (directory / 'queries.fna').read_bytes()
    query_hash = hashlib.sha256(raw).hexdigest()
    queries = [line[1:].split()[0] for line in raw.decode().splitlines() if line.startswith('>')]
    assert len(set(queries)) == len(queries)
    protocols = {}
    hashes = {}
    for method in (left, right):
        d = directory / method
        complete = json.loads((d / 'complete.json').read_text())
        resource = json.loads((d / 'resources.json').read_text())
        protocol = json.loads((d / 'protocol.json').read_text())
        assert resource['status'] == 'complete' and complete['queries'] == len(queries)
        assert protocol['input_manifest']['query_sha256'] == query_hash
        digest = hashlib.sha256((d / 'alignments.tsv').read_bytes()).hexdigest()
        assert digest == complete['output_sha256']
        protocols[method] = protocol
        hashes[method] = digest
    assert protocols[left]['input_manifest'] == protocols[right]['input_manifest']
    result, per_query, only = compare(directory / left / 'alignments.tsv', directory / right / 'alignments.tsv', queries)
    result.update(left_method=left, right_method=right, query_sha256=query_hash,
                  input_hashes=hashes, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  search_source_commits={name: x['source_commit'] for name, x in protocols.items()},
                  commands={name: x['command'] for name, x in protocols.items()},
                  limits=['Comparison preserves all 19 original output fields and duplicate multiplicity; row order is ignored.',
                          'Every frozen query remains in the denominator, including queries without returned alignments.',
                          'Raw-output equality proves parity on this finite panel only; it is not biological ground truth.',
                          'Runtime comparison uses the same worker and reused reference cache, not randomized cold-cache trials.'])
    output.mkdir(exist_ok=False)
    with (output / 'per-query.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(per_query[0]), delimiter='\t')
        w.writeheader(); w.writerows(per_query)
    with gzip.open(output / 'raw-differences.tsv.gz', 'wt') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['side', 'count'] + FIELDS)
        for label, rows in zip((left, right), only):
            for row, count in sorted(rows.items()):
                w.writerow([label, count, *row])
    (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('commands', 'limits')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('directory', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--left-method', default='tblastx-word4')
    p.add_argument('--right-method', default='tblastx-word4-one-batch')
    a = p.parse_args()
    main(a.directory, a.left_method, a.right_method, a.output)

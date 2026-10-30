"""Describe nested depths using accepted aggregate counts, without new calls."""
import argparse
import csv
import json
from pathlib import Path

METRICS = ('production_rows', 'protein_rows', 'production_only', 'protein_only',
           'production_joint_exemplars', 'protein_variant_joint_exemplars')

def compare(rows, depth_runs):
    index = {(r['sample'], int(r['depth_million_pairs'])): r for r in rows}
    assert len(index) == len(rows), 'Duplicate run/depth rows'
    result = []
    for sample in depth_runs:
        shallow, deep = index[(sample, 5)], index[(sample, 25)]
        row = dict(sample=sample, shallow_status=shallow['status'], deep_status=deep['status'])
        row['status'] = 'complete' if shallow['status'] == deep['status'] == 'complete' else 'pending'
        if 'blocked' in (shallow['status'], deep['status']):
            row['status'] = 'blocked'
        if row['status'] == 'complete':
            for metric in METRICS:
                lo, hi = int(shallow[metric]), int(deep[metric])
                assert lo >= 0 and hi >= 0
                row[metric+'_5m'] = lo
                row[metric+'_25m'] = hi
                row[metric+'_5m_per_million_input_pairs'] = lo/5
                row[metric+'_25m_per_million_input_pairs'] = hi/25
                row[metric+'_count_ratio'] = hi/lo if lo else None
        result.append(row)
    return result

def main():
    p = argparse.ArgumentParser()
    for key in ('aggregate', 'selection', 'out'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    with a.aggregate.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    result = compare(rows, json.loads(a.selection.read_text())['depth_runs'])
    a.out.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(k for row in result for k in row))
    with (a.out/'nested-depth.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader(); writer.writerows(result)
    (a.out/'nested-depth.json').write_text(json.dumps(dict(
        planned_pairs=len(result), complete_pairs=sum(r['status']=='complete' for r in result),
        interpretation=[
            'Each row compares nested input subsets of one sequencing run, not independent observations.',
            'Counts are pipeline output records or production duplicate exemplars, not confirmed biological detections.',
            'A difference of totals is not the number of newly observed records or distinct origins.',
            'Five times as many input pairs does not imply five times as many output records or exemplars.',
            'A zero shallow count has no defined count ratio; blank is not zero.',
            'Pending or blocked comparisons are not included as zero results.',
            'Timing scaling requires the separate runtime-context table because cache state and concurrent work differ.'
        ]), indent=2)+'\n')

if __name__ == '__main__':
    main()

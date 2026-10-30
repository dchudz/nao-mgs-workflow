"""Summarize existing acquisition receipts once per source run, preserving unknowns."""
import argparse, collections, csv, datetime, hashlib, json, statistics
from pathlib import Path

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    inventory = args.study / 'results/aggregate/runs.tsv'
    with inventory.open() as handle:
        accepted = [r for r in csv.DictReader(handle, delimiter='\t') if r['status'] == 'complete']
    grouped = collections.defaultdict(list)
    inputs = {str(inventory): sha(inventory)}
    for row in accepted:
        path = args.study / 'results' / row['label'] / 'input-acquisition.json'
        receipt = json.loads(path.read_text())
        inputs[str(path)] = sha(path)
        grouped[receipt['sample']].append(receipt)
    per_run = []
    keys = ['source_sha256', 'population_pairs', 'source_bytes', 'seconds_excluding_initial_download',
            'vdb_validation_seconds', 'download_seconds']
    for sample, receipts in sorted(grouped.items()):
        first = receipts[0]
        assert all(all(r.get(k) == first.get(k) for k in keys) for r in receipts), sample
        elapsed = first['seconds_excluding_initial_download']
        validation = first.get('vdb_validation_seconds')
        assert elapsed >= 0 and (validation is None or 0 <= validation <= elapsed)
        download = first.get('download_seconds')
        assert download is None or download >= 0
        per_run.append(dict(sample=sample, accepted_comparisons=len(receipts), source_bytes=first['source_bytes'],
            source_pairs=first['population_pairs'], recorded_download_seconds=download,
            final_acquisition_excluding_download_seconds=elapsed,
            validation_seconds_within_acquisition=validation,
            remaining_acquisition_seconds=None if validation is None else elapsed-validation))
    metrics = []
    for field in ['recorded_download_seconds', 'final_acquisition_excluding_download_seconds',
                  'validation_seconds_within_acquisition', 'remaining_acquisition_seconds']:
        values = [r[field] for r in per_run if r[field] is not None]
        metrics.append(dict(metric=field, source_runs=len(per_run), recorded_runs=len(values),
            missing_runs=len(per_run)-len(values), recorded_hours=sum(values)/3600,
            median_recorded_minutes=statistics.median(values)/60 if values else None))
    args.out.mkdir(parents=True, exist_ok=True)
    for name, records in [('acquisition-per-run', per_run), ('acquisition-summary', metrics)]:
        with (args.out / (name+'.tsv')).open('w') as handle:
            writer=csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
            writer.writeheader(); writer.writerows(records)
    proof = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        accepted_comparisons=len(accepted), unique_source_runs=len(per_run),
        analysis_sha256=sha(Path(__file__)), input_hashes=inputs,
        outputs={name:sha(args.out/name) for name in ['acquisition-per-run.tsv','acquisition-summary.tsv']},
        definitions={'unit':'One acquisition per source run, shared across arms and nested depths.',
        'unknown':'Missing download timing remains blank and is excluded from recorded-time totals.',
        'overlap':'Validation is included in final acquisition time; do not add it again.',
        'limits':'Successful final acquisition plus documented download receipts is not all failed/interrupted work or worker billing. The worker ledger includes those resource intervals.'})
    (args.out/'acquisition-provenance.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({'accepted_comparisons':len(accepted), 'source_runs':len(per_run), 'metrics':metrics}))

if __name__ == '__main__':
    main()

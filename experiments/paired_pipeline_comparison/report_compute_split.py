"""Summarize measured stage durations separately from worker cost and acquisition."""
import argparse
import csv
import json
import statistics
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, required=True)
    study = parser.parse_args().study
    result = study/'results'
    with (result/'aggregate/runs.tsv').open() as stream:
        completed = [r for r in csv.DictReader(stream, delimiter='\t') if r['status']=='complete']
    groups = {}
    for unit in completed:
        with (result/unit['label']/'comparison-report/stage-wall-times.tsv').open() as stream:
            stages = list(csv.DictReader(stream, delimiter='\t'))
        measured = sum(float(r['wall_seconds']) for r in stages)
        assert abs(measured/60-float(unit['stage_minutes'])) < 1e-6
        assert len(stages) == len({r['stage'] for r in stages})
        for row in stages:
            groups.setdefault((int(unit['depth_million_pairs']),row['stage']),[]).append(float(row['wall_seconds'])/60)
    rows = [dict(depth_million_pairs=depth,stage=stage,units=len(values),median_minutes=statistics.median(values),minimum_minutes=min(values),maximum_minutes=max(values),total_minutes=sum(values)) for (depth,stage),values in sorted(groups.items())]
    output = result/'aggregate/compute-split.tsv'
    if rows:
        with output.open('w') as stream:
            writer=csv.DictWriter(stream, fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    text = ['# Measured stage time', '', 'These are completed-unit wall times, not CPU times or a bill. Cold reference work, validation order, and sampled records can affect them. The 5M and 25M strata are nested and are kept separate.', '', '| Pairs | Stage | Units | Median minutes | Range minutes |', '|---|---|---:|---:|---:|']
    for row in rows:
        text.append(f"| {row['depth_million_pairs']}M | {row['stage']} | {row['units']} | {row['median_minutes']:.2f} | {row['minimum_minutes']:.2f}–{row['maximum_minutes']:.2f} |")
    text.extend(['', 'Acquisition of the entire public source archive, input sampling, boot, preservation, retries and idle periods are outside this table. Acquisition is shared by both arms and by nested depths. The execution-overhead receipts account for wrapper time outside the listed stages; they do not attribute that residual to a search.', '', 'The baseline run also produces QC/profiling outputs reused by the variant. Adding the protein stages to baseline time is the incremental cost of this paired experiment, not a standalone deployment cost ratio. A standalone variant would still need common QC/profiling work, whose cost is not isolated by these top-level stages.', '', 'Use runtime-context.tsv for measured overlaps with reference-storage checks and the disjoint worker-interval ledger for compute-cost bounds. No overlap duration is subtracted here.', ''])
    (result/'aggregate/compute-split.md').write_text('\n'.join(text))
    print(json.dumps({'completed_units':len(completed),'stage_groups':len(rows)}))


if __name__=='__main__':main()

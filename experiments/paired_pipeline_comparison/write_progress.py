"""Render a short reading guide from accepted reports; do not classify reads."""
import argparse
import csv
import datetime
import json
from pathlib import Path

def read(path):
    return json.loads(path.read_text())

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--study', type=Path, required=True)
    a = p.parse_args(); w = a.study
    with (w/'results/aggregate/runs.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    complete = [r for r in rows if r['status']=='complete']
    pilot = [r for r in rows if r['feasibility_run']=='True']
    pilot_done = [r for r in pilot if r['status']=='complete']
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    expansion = w/'runtime-provenance/expansion/accepted-expansion-plan.json'
    expanded = expansion.exists() and read(expansion)['status']=='accepted-for-execution'
    execution = (f'{len(rows)-len(complete)} of the full cohort comparisons remain uncollected. The feasibility gate passed and the fixed expansion was accepted within the budget. ' if expanded else f'{len(pilot)-len(pilot_done)} feasibility comparisons remain uncollected. Expansion is conditional on completed outputs and measured costs. ')
    execution += 'No sample is chosen or discarded because of its biological findings. Uncollected does not mean a zero result.'
    status_path = w/'execution-status.json'
    halted = status_path.exists() and read(status_path)['status']=='halted'
    if halted:
        execution = (f'{len(rows)-len(complete)} planned comparisons are unfinished. Execution is halted; no pipeline jobs remain running. '
                     'This was an assistant scope decision, not budget exhaustion or a failed correctness gate. '
                     'The completed subset is not a random completion sample and does not establish population-level performance. '
                     'Unfinished units are not zero results. See [the shutdown and preservation record](../runtime-provenance/preservation-halt/README.md).')
    lines = ['# Frozen paired pipeline comparison', '', f'Report snapshot: {now}.', '',
        f'{len(complete)} of {len(rows)} planned comparisons are collected and accepted, including {len(pilot_done)} of {len(pilot)} feasibility comparisons. '
        f'The full frozen design contains {len(rows)} comparisons: 50 at 5 million pairs and 10 nested comparisons at 25 million pairs. '
        + ('Execution is halted; unfinished units are retained as incomplete.' if halted else 'A missing local report can still be running remotely.'), '',
        'Both arms use the same selected input pairs, source reference collection, production duplicate handling, '
        'validation mode and validation cap. The variant changes the upstream screen and alignment representation. '
        'This does not isolate alignment as the only cause of differences.', '',
        '## Completed output inventory', '',
        'These are final output records, not confirmed detections or species calls. '
        '“Only” compares pair identifiers between arms at the same depth. It does not compare shallow and deep identifiers. Reassigned shared records occur in both arms with different labels; these are separate from added or lost records.', '',
        '| Run | Input pairs (million) | Production records | Protein records | Production only | Protein only | Shared records reassigned | Stage minutes | Driver minutes |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in complete:
        overhead_path = w/'results'/r['label']/'execution-overhead.json'
        overhead = read(overhead_path) if overhead_path.exists() else {}
        receipt = overhead.get('pair_receipt')
        minutes = f"{receipt['seconds']/60:.1f}" if receipt else 'unrecorded'
        lines.append(f"| {r['sample']} | {r['depth_million_pairs']} | {int(r['production_rows']):,} | "
                     f"{int(r['protein_rows']):,} | {int(r['production_only']):,} | {int(r['protein_only']):,} | "
                     f"{int(r['reassigned']):,} | {float(r['stage_minutes']):.1f} | {minutes} |")
    lines += ['', 'Stage minutes omit preflight, inter-stage backups and driver overhead. Driver minutes include those, '
        'but omit acquisition and the final verified archive. Worker uptime captures all of them. '
        'Early replica timings overlap storage initialization; see [runtime context](aggregate/runtime-context.tsv). '
        'The same input QC/profiling work is shared between arms, so these are not standalone deployment timings.', '',
        '## Interpretation', '',
        '- The known-origin artificial control passed 24 positive pairs, one negative pair, and 12 expected duplicate groups using the production duplicate container.',
        '- Every accepted unit has both arms complete, equal downstream parameters, output checks, and a checksum-verified private S3 archive.',
        '- Bracken abundance is unavailable because the retained taxonomy/report combination has no requested D-rank rows. These are accepted assignment comparisons, not complete abundance reports.',
        '- Clade counts precede validation. Validation annotates sampled records and does not turn all output records into confirmed positives.',
        '- Standard validation agreement allows a parent label to agree with a more specific child label. It is not exact species agreement or accuracy.',
        '- Output additions, losses and changes have no complete environmental truth labels. More output is not established to be better.',
        '- Equal downstream settings do not imply matched upstream scoring thresholds or false-positive operating points.',
        '- Nested depth comparisons reuse observations from the same run. Differences of exemplar totals are not counts of new independent origins.', '',
        'See [the interpretation record](../comparison-interpretation.md), [nested depth table](aggregate/nested-depth.tsv), '
        '[measured stage split](aggregate/compute-split.md), and [full acceptance inventory](aggregate/runs.tsv). Individual run folders retain the standard comparison reports, source commit, '
        'commands, stage timings, acquisition receipts and private archive checksum.', '',
        '## Remaining execution', '',
        execution, '',
        'The existing budget is $500. The working allocation reserves $200 for new compute; '
        'posted historical charges, billing lag and noncompute retention are accounted for separately in [budget.json](../budget.json). '
        'Observed resource-time ledgers live in `../runtime-provenance/`. They are estimates and ceilings, not invoice reconciliation.', '']
    (w/'results/progress.md').write_text('\n'.join(lines))

if __name__ == '__main__':
    main()

"""Measure query routing against completed comprehensive searches; never infer runtime savings."""
import argparse, collections, csv, hashlib, json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('analysis', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--left', default='blastn')
p.add_argument('--right', default='tblastx-word4')
a = p.parse_args()
source = a.analysis / 'taxonomy-description/per-mate.tsv'
rows = list(csv.DictReader(source.open(), delimiter='\t'))
a.output.mkdir(parents=True, exist_ok=False)
rules = {
    'no_alignment': lambda x: x['natural_status'] == 'no_alignment',
    'not_viral_only': lambda x: x['natural_status'] != 'viral_only',
    'not_viral_species_rank': lambda x: x['natural_status'] != 'viral_only' or x['pipeline_natural_lca_rank'] != 'species',
    'all_queries': lambda x: True,
}
definitions = {
    'natural_verdict': ('natural_status', 'pipeline_natural_lca'),
    'combined_verdict': ('status', 'pipeline_combined_lca'),
    'near_best_taxon_set': ('near_best_taxids',),
}
details = []
summary = []
for threshold in ('default_E10', 'strict_E1e-5'):
    left = {x['query']: x for x in rows if x['method'] == a.left and x['threshold'] == threshold}
    right = {x['query']: x for x in rows if x['method'] == a.right and x['threshold'] == threshold}
    assert left and left.keys() == right.keys()
    for rule, select in rules.items():
        routed = {q for q, x in left.items() if select(x)}
        for label, fields in definitions.items():
            changed = {q for q in left if tuple(left[q][f] for f in fields) != tuple(right[q][f] for f in fields)}
            captured = changed & routed
            missed = changed - routed
            assert len(captured) + len(missed) == len(changed)
            if rule == 'all_queries': assert not missed
            summary.append(dict(threshold=threshold, rule=rule, comparison=label, queries=len(left), routed_queries=len(routed), unrouted_queries=len(left)-len(routed), changed_queries=len(changed), captured_changes=len(captured), missed_changes=len(missed), missed_pairs=len({left[q]['pair'] for q in missed})))
            for q in sorted(changed):
                details.append(dict(threshold=threshold, rule=rule, comparison=label, query=q, pair=left[q]['pair'], stratum=left[q]['stratum'], routed=q in routed, left_status=left[q]['natural_status'], right_status=right[q]['natural_status'], left_lca=left[q]['pipeline_natural_lca'], right_lca=right[q]['pipeline_natural_lca'], left_best_evalue=left[q]['best_evalue'], right_best_evalue=right[q]['best_evalue']))
for name, data in [('routing-counts.tsv', summary), ('changed-read-routing.tsv', details)]:
    with (a.output / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]) if data else ['none'], delimiter='\t'); w.writeheader(); w.writerows(data)
result = dict(left=a.left, right=a.right, input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), counts=summary, limits=[
    'Retrospective routing analysis of completed full-reference results; no additional search and no query verdict is transferred between reads.',
    'Rules were selected after seeing the small pilot. They are descriptive comparisons, not independently tuned or validated policies.',
    'Counts describe disagreements, not biological truth, sensitivity, false positives, or validated accuracy.',
    'Query-count reduction does not imply proportional CPU, elapsed-time or cost savings; reference scans and batching must be measured.',
    'The species-rank rule tests the reported LCA rank literally; it also routes any below-species rank and is not a species-resolution classifier.',
    'Natural verdict excludes artificial constructs as in the source annotation. Combined verdict and near-best taxon sets are reported separately.',
    'Every query remains in the simulated cascade: routed queries take the translated verdict and all other queries retain the nucleotide verdict.',
])
(a.output / 'summary.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps([x for x in summary if x['threshold']=='strict_E1e-5' and x['comparison']=='natural_verdict'], indent=2))

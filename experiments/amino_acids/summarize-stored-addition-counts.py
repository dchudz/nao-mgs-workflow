"""Aggregate saved flags and annotations without reading or reprocessing sequences."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('heldout', type=Path)
parser.add_argument('annotation', type=Path)
parser.add_argument('output', type=Path)
a = parser.parse_args()
paths = {name: a.heldout / name for name in ['added-pairs.tsv', 'sample-metrics.tsv', 'summary.json']}
paths['annotation'] = a.annotation
with paths['added-pairs.tsv'].open() as f:
    added = [r for r in csv.DictReader(f, delimiter='\t') if r['arm'] == 'union']
with paths['sample-metrics.tsv'].open() as f:
    metrics = [r for r in csv.DictReader(f, delimiter='\t') if r['arm'] == 'union']
summary = json.loads(paths['summary.json'].read_text())
with gzip.open(a.annotation, 'rt') as f:
    annotation = {r['taxid']: r for r in csv.DictReader(f, delimiter='\t')}
assert len({(r['sample'], r['seq_id']) for r in added}) == len(added)
assert len(metrics) == len({r['sample'] for r in metrics}) == summary['completed_samples']
by_sample = collections.Counter(r['sample'] for r in added)
assert set(by_sample) <= {r['sample'] for r in metrics}
assert all(by_sample[r['sample']] == int(r['added_pairs']) for r in metrics)
assert len(added) == summary['arms']['union']['totals']['added_pairs']
assert all(r['joint_exemplar'] in ['True', 'False'] for r in added)
values = sorted(by_sample[r['sample']] for r in metrics)
by_species = collections.Counter()
unresolved = collections.Counter()
for r in added:
    node = annotation.get(r['taxid'])
    species = (node['taxid'] if node['rank'] == 'species' else node['taxid_species']) if node else None
    if species in annotation and annotation[species]['rank'] == 'species':
        by_species[species] += 1
    else:
        unresolved[node['rank'] if node else 'missing_taxid'] += 1
assert sum(by_species.values()) + sum(unresolved.values()) == len(added)
result = {
    'scope': 'Local aggregation of already saved per-pair flags, sample counts, and frozen taxonomy annotations. No sequence data, new search, duplicate rerun, or species calling.',
    'source_files': {k: {'path': str(p.resolve()), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for k, p in paths.items()},
    'added_pairs': len(added),
    'added_pairs_marked_joint_exemplars_in_full_union': sum(r['joint_exemplar'] == 'True' for r in added),
    'other_added_pairs': sum(r['joint_exemplar'] == 'False' for r in added),
    'joint_exemplar_definition': 'Saved flag from analyze-cohort.py: seq_id == prim_align_dup_exemplar == sim_dup_exemplar, within the original full-union run for each sample.',
    'distinct_duplicate_groups_touched_by_additions': None,
    'exemplar_limit': 'The saved table has a boolean flag, not representative IDs. It does not give the number of all duplicate groups touched by additions, an additions-only reclustering result, or distinct loci across samples. Representative selection can change when records are removed.',
    'samples': {'cohort': len(metrics), 'with_additions': len(by_sample), 'without_additions': len(metrics)-len(by_sample), 'median_all_samples': statistics.median(values), 'mean_all_samples': sum(values)/len(values), 'maximum': max(values), 'top_10_pairs': sum(values[-10:]), 'top_10_share': sum(values[-10:])/len(added), 'top_100_pairs': sum(values[-100:]), 'top_100_share': sum(values[-100:])/len(added), 'nearest_rank_percentiles': {str(q): values[math.ceil(q*len(values))-1] for q in [.9,.95,.99]}},
    'species': {'distinct_resolved_species': len(by_species), 'resolved_pairs': sum(by_species.values()), 'unresolved_pairs': sum(unresolved.values()), 'unresolved_by_rank': dict(unresolved), 'largest_species_pairs': max(by_species.values()), 'largest_species_share_all_additions': max(by_species.values())/len(added), 'largest_species_share_resolved': max(by_species.values())/sum(by_species.values()), 'definition': 'Use the frozen annotation taxid_species field, or the assigned taxid itself when its rank is species. Assignments above species remain unresolved; no evidence is redistributed.'},
    'interpretation': 'The overlap is exact descriptive agreement on this cohort, not a precision estimate from independent events. Per-sample exemplars can recur at the same loci across samples. Sample, location, and species concentration remain relevant.'
}
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k in ['added_pairs', 'added_pairs_marked_joint_exemplars_in_full_union', 'samples', 'species']}))

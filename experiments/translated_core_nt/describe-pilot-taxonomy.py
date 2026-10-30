"""Separate artificial constructs from natural competitors without changing frozen calls."""
import argparse
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--directory', type=Path, required=True)
p.add_argument('--analysis', default='analysis')
a = p.parse_args()
r = Path('/home/ec2-user/experiment')
source = a.directory / a.analysis
out = source / 'taxonomy-description'
out.mkdir(exist_ok=False)
parents, ranks = {}, {}
for line in (r / 'reference/results/taxonomy-nodes.dmp').open():
    f = [x.strip() for x in line.split('|')]
    parents[int(f[0])], ranks[int(f[0])] = int(f[1]), f[2]
cache = {}
def lineage(t):
    if t not in cache:
        original = t
        found = set()
        while t in parents and t not in found:
            found.add(t)
            t = parents[t]
        cache[original] = found
        return found
    return cache[t]
def category(t):
    path = lineage(t)
    if not path:
        return 'unknown'
    for ancestor, name in [(81077, 'artificial'), (10239, 'viral'), (2, 'bacterial'), (2157, 'archaeal'), (2759, 'eukaryotic')]:
        if ancestor in path:
            return name
    return 'other_known'
rows = list(csv.DictReader((source / 'per-mate.tsv').open(), delimiter='\t'))
wanted = set()
for x in rows:
    taxa = {int(t) for t in x['near_best_taxids'].split(';') if t}
    groups = collections.defaultdict(set)
    for t in taxa:
        groups[category(t)].add(t)
    natural = taxa - groups['artificial']
    if not natural:
        status = 'artificial_only' if taxa else 'no_alignment'
    elif groups['unknown']:
        status = 'unknown_taxonomy'
    else:
        viral = bool(groups['viral'])
        nonviral = bool(natural - groups['viral'])
        status = 'mixed_viral_natural_nonviral' if viral and nonviral else 'viral_only' if viral else 'natural_nonviral_only'
    x.update(natural_status=status, artificial_taxids=';'.join(map(str, sorted(groups['artificial']))), natural_nonviral_taxids=';'.join(map(str, sorted(natural - groups['viral']))))
    wanted.update(taxa)
    if x['pipeline_natural_lca'] != 'NA':
        wanted.add(int(x['pipeline_natural_lca']))
names = {}
for line in (r / 'reference/results/taxonomy-names.dmp').open():
    f = [x.strip() for x in line.split('|')]
    if f[3] == 'scientific name' and int(f[0]) in wanted:
        names[int(f[0])] = f[1]
for x in rows:
    t = int(x['pipeline_natural_lca']) if x['pipeline_natural_lca'] != 'NA' else None
    x.update(pipeline_natural_lca_name=names.get(t, ''), pipeline_natural_lca_rank=ranks.get(t, ''))
methods = sorted({x['method'] for x in rows})
thresholds = sorted({x['threshold'] for x in rows})
lookup = {(x['method'], x['threshold'], x['query']): x for x in rows}
queries = sorted({x['query'] for x in rows})
changes, comparisons = [], {}
for threshold in thresholds:
    for left, right in itertools.combinations(methods, 2):
        stats = collections.Counter()
        for q in queries:
            l, rr = lookup[left, threshold, q], lookup[right, threshold, q]
            before, after = l['pipeline_natural_lca'], rr['pipeline_natural_lca']
            relation = 'unchanged' if before == after else 'gained_alignment' if before == 'NA' else 'lost_alignment' if after == 'NA' else 'broader_LCA' if int(after) in lineage(int(before)) else 'more_specific_LCA' if int(before) in lineage(int(after)) else 'different_lineage'
            stats[relation] += 1
            stats['natural_status_changed'] += int(l['natural_status'] != rr['natural_status'])
            if relation != 'unchanged' or l['natural_status'] != rr['natural_status']:
                changes.append(dict(query=q, stratum=l['stratum'], threshold=threshold, left=left, right=right, left_status=l['natural_status'], right_status=rr['natural_status'], left_lca=before, right_lca=after, left_name=l['pipeline_natural_lca_name'], right_name=rr['pipeline_natural_lca_name'], relation=relation))
        comparisons[left + ' -> ' + right + '/' + threshold] = dict(stats)
raw_caps = {}
for method in methods:
    subjects = collections.defaultdict(set)
    with (a.directory / method / 'alignments.tsv').open() as f:
        for values in csv.reader(f, delimiter='\t'):
            subjects[values[0]].add(values[1])
    raw_caps[method] = dict(queries_with_at_least250_returned_subjects=sum(len(v) >= 250 for v in subjects.values()), maximum_returned_subjects=max(map(len, subjects.values()), default=0))
for name, data in [('per-mate.tsv', rows), ('natural-LCA-changes.tsv', changes), ('taxa.tsv', [dict(taxid=t, name=names.get(t, ''), rank=ranks.get(t, ''), category=category(t)) for t in sorted(wanted)])]:
    with (out / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]) if data else ['none'], delimiter='\t')
        w.writeheader()
        w.writerows(data)
summary = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), original_per_mate_sha256=hashlib.sha256((source / 'per-mate.tsv').read_bytes()).hexdigest(), natural_status_counts={method + '/' + threshold: dict(collections.Counter(x['natural_status'] for x in rows if x['method'] == method and x['threshold'] == threshold)) for method in methods for threshold in thresholds}, rows_with_artificial_hits={method + '/' + threshold: sum(bool(x['artificial_taxids']) for x in rows if x['method'] == method and x['threshold'] == threshold) for method in methods for threshold in thresholds}, natural_LCA_comparisons=comparisons, subject_cap_diagnostic=raw_caps, limits=['Supplement annotates existing outputs; search results, score filtering, and production LCA calculations are unchanged.', 'Original mixed_viral_nonviral status includes artificial constructs. It must not be interpreted as natural bacterial/human competition without this distinction.', 'Taxonomic comparisons describe returned, filtered hits. The250subject reporting cap can omit competing references.', 'Broader LCA is taxonomic resolution loss, not a demonstrated false viral detection. Different methods and missing translated hits do not by themselves establish biological truth.'])
(out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary))

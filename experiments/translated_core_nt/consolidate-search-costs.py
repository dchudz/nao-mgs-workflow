"""Tabulate only completed measured searches, retaining unequal-batch caveats."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('root', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
helper = Path(__file__).with_name('compare-search-costs.py')
spec = importlib.util.spec_from_file_location('costs', helper)
costs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(costs)
rows, pending, inputs = [], [], {}


def read(path):
    data = path.read_bytes()
    inputs[str(path)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)


for panel, methods in [('pilot24', ['megablast', 'blastn', 'tblastx', 'query-tblastn-fast', 'tblastx-word4']),
                       ('protein232', ['megablast', 'blastn', 'tblastx-word4'])]:
    for method in methods:
        d = a.root / panel / method
        if not (d / 'complete.json').exists():
            pending.append(panel + '/' + method)
            continue
        done = read(d / 'complete.json')
        resource = read(d / 'resources.json')
        protocol = read(d / 'protocol.json')
        assert resource['status'] == 'complete'
        count = protocol['input_manifest']['queries']
        timing = costs.timing(d / 'time.txt')
        inputs[str(d / 'time.txt')] = hashlib.sha256((d / 'time.txt').read_bytes()).hexdigest()
        rows.append(dict(panel=panel, method=method, reference='core_nt', query_mates=count,
                         query_sha256=protocol['input_manifest']['query_sha256'],
                         wall_seconds=timing['search_wall_seconds'], CPU_seconds=timing['CPU_seconds'],
                         peak_RSS_KiB=timing['max_RSS_KiB'], threads=16, instance_type='x2iedn.4xlarge',
                         allocated_USD=timing['search_wall_seconds'] / 3600 * 3.6365,
                         command=protocol['command'], source_commit=protocol['source_commit']))

context = read(a.root / 'seed-nr-standalone/context-cost-and-mode.json')
original = read(a.root / 'completed-search-cost-comparison.json')
pooled_command = next(x['command'] for x in original['searches'] if x['method'] == 'DIAMOND-sensitive')
standalone_protocol = read(a.root / 'seed-nr-standalone/search-protocol.json')
for name, x in context['searches'].items():
    rows.append(dict(panel=name, method='DIAMOND-sensitive/' + ','.join(x['algorithms']), reference='NR',
                     query_mates=x['queries'], query_sha256=x['query_sha256'],
                     wall_seconds=x['GNU_wall_seconds'], CPU_seconds=x['GNU_user_seconds'] + x['GNU_system_seconds'],
                     peak_RSS_KiB=x['GNU_max_RSS_KiB'], threads=32, instance_type='r7a.8xlarge',
                     allocated_USD=x['allocated_compute_EBS_IPv4_USD'],
                     command=pooled_command if name == 'combined-validation' else standalone_protocol['command']))
fixed = a.root / 'seed-nr-forced-double-index-trace'
if (fixed / 'complete.json').exists():
    read(fixed / 'complete.json')
    assert read(fixed / 'resources.json')['status'] == 'complete'
    protocol = read(fixed / 'search-protocol.json')
    log = fixed / 'nr-search.log'
    inputs[str(log)] = hashlib.sha256(log.read_bytes()).hexdigest()
    assert 'Algorithm: Double-indexed' in log.read_text()
    timing = costs.timing(log)
    rate = 2.4344 + (1536 * .08 + 3000 * .005 + 875 * .04) / 720 + .005
    rows.append(dict(panel=fixed.name, method='DIAMOND-sensitive/forced-Double-indexed', reference='NR',
                     query_mates=protocol['queries'], query_sha256=protocol['query_sha256'],
                     wall_seconds=timing['search_wall_seconds'], CPU_seconds=timing['CPU_seconds'],
                     peak_RSS_KiB=timing['max_RSS_KiB'], threads=32, instance_type='r7a.8xlarge',
                     allocated_USD=timing['search_wall_seconds'] / 3600 * rate,
                     command=protocol['command'], source_commit=protocol['source_commit']))
else:
    pending.append(fixed.name)
for row in rows:
    row['allocated_USD_per_query'] = row['allocated_USD'] / row['query_mates']
    row['CPU_seconds_per_query'] = row['CPU_seconds'] / row['query_mates']
    row['reference_file_bytes'] = 318551974004 if row['reference'] == 'core_nt' else 918566558190
    row['new_persistent_search_index_built'] = False

lookup = {(x['panel'], x['method']): x for x in rows}
ratios = []
for method in ['tblastx', 'query-tblastn-fast', 'tblastx-word4']:
    left = lookup['pilot24', method]
    for right in rows:
        if right['reference'] == 'NR':
            ratios.append(dict(numerator='pilot24/' + method, denominator=right['panel'] + '/' + right['method'],
                               observed_per_query_allocation_ratio=left['allocated_USD_per_query'] / right['allocated_USD_per_query']))
result = dict(searches=rows, searches_without_local_completion_receipt=pending,
              unequal_experiment_cost_ratios=ratios, input_sha256=inputs,
              source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              limits=['All figures are completed-search wall-time allocations, already included in lifetime spend.',
                      'Different query pools, references, hardware, cache states, scoring and automatic modes prevent an isolated algorithm-speed interpretation.',
                      'The 43-mate and 464-mate core_nt panels are not nested and have different selected strata.',
                      'Full reference inventories are reused native BLAST databases. No persistent search index was newly built; reference preparation and transfer are outside these search timers.',
                      'No peak temporary-disk measurement was collected. Reference-file bytes are not total worker disk allocation.',
                      'Peak RSS includes mapped reference pages and is not a measured minimum RAM requirement.',
                      'Pooled NR time is for all 817744 queries, not an independently timed 200-sample-only run.',
                      'Incomplete searches and the intentionally interrupted control are excluded from this completed-search table but charged in the lifetime ledger.'])
a.output.mkdir(exist_ok=False)
(a.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
columns = ['panel', 'method', 'reference', 'query_mates', 'wall_seconds', 'CPU_seconds', 'threads',
           'instance_type', 'peak_RSS_KiB', 'allocated_USD', 'allocated_USD_per_query', 'reference_file_bytes']
with (a.output / 'measurements.tsv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    w.writerows(rows)
print(json.dumps({'completed_searches': len(rows), 'ratios': ratios}))

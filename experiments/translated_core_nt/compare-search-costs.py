"""Compare completed search allocations with explicit batching/reference limitations."""
import argparse
import datetime
import hashlib
import json
import re
from pathlib import Path


def timing(path):
    text = Path(path).read_text()
    assert re.search(r'Exit status:\s*0\s*$', text), 'Successful GNU time record required'
    labels = {'user_seconds': 'User time (seconds)', 'system_seconds': 'System time (seconds)', 'max_RSS_KiB': 'Maximum resident set size (kbytes)'}
    result = {}
    for key, label in labels.items():
        match = re.search(re.escape(label) + r':\s*([\d.]+)', text)
        assert match, label
        result[key] = float(match[1])
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', text)
    assert match
    result['search_wall_seconds'] = 0
    for part in match[1].split(':'):
        result['search_wall_seconds'] = result['search_wall_seconds'] * 60 + float(part)
    result['CPU_seconds'] = result['user_seconds'] + result['system_seconds']
    return result


def row(name, reference, queries, rate, metrics, command, query_hash):
    cost = metrics['search_wall_seconds'] / 3600 * rate
    return dict(method=name, reference=reference, query_mates=queries, query_sha256=query_hash, command=command, modeled_hourly_allocation_USD=rate, search_allocation_USD=cost, search_allocation_USD_per_query=cost / queries, CPU_seconds_per_query=metrics['CPU_seconds'] / queries, **metrics)


def main(a):
    pilot, nr = Path(a.pilot), Path(a.nr)
    rates = {'primary': 3.3345 + (2048 * .08 + 3000 * .005 + 875 * .04) / 720 + .005, 'NR': 2.4344 + (1280 * .08 + 3000 * .005 + 875 * .04) / 720 + .005}
    rows = []
    for method in ['megablast', 'blastn', 'tblastx']:
        d = pilot / method
        complete = json.loads((d / 'complete.json').read_text())
        protocol = json.loads((d / 'protocol.json').read_text())
        resources = json.loads((d / 'resources.json').read_text())
        assert resources['status'] == 'complete'
        result = row(method, 'full frozen core_nt nucleotide reference', complete['queries'], rates['primary'], timing(d / 'time.txt'), protocol['command'], protocol['input_manifest']['query_sha256'])
        result.update(runner_wall_seconds=resources['wall_seconds'], source_commit=protocol['source_commit'])
        rows.append(result)
    raw_receipt = json.loads((nr / 'search-raw-complete.json').read_text())
    for name in ['nr-search.log', 'search-protocol.json']:
        assert (nr / name).stat().st_size == raw_receipt['files_bytes'][name]
    assert raw_receipt['files_bytes']['nr.tsv'] > 0
    protocol = json.loads((nr / 'search-protocol.json').read_text())
    rows.append(row('DIAMOND-sensitive', 'full frozen NR protein reference', protocol['queries'], rates['NR'], timing(nr / 'nr-search.log'), protocol['command'], protocol['query_sha256']))
    # Preserve the pre-expansion allocation, then charge added capacity only
    # for the overlap with the completed search's recorded interval.
    expansion = json.loads((nr / 'nr-volume-expansion-20261022.json').read_text())
    expanded_at = datetime.datetime.fromisoformat(expansion['modification']['VolumeModification']['StartTime'])
    search_started = datetime.datetime.fromisoformat(protocol['started_utc'])
    search_ended = search_started + datetime.timedelta(seconds=rows[-1]['search_wall_seconds'])
    overlap = max(0, (search_ended - max(expanded_at, search_started)).total_seconds())
    extra = overlap / 3600 * expansion['additional_storage_USD_per_hour']
    rows[-1]['base_hourly_allocation_USD'] = rates['NR']
    rows[-1]['added_storage_overlap_seconds'] = overlap
    rows[-1]['added_storage_allocation_USD'] = extra
    rows[-1]['search_allocation_USD'] += extra
    rows[-1]['search_allocation_USD_per_query'] = rows[-1]['search_allocation_USD'] / rows[-1]['query_mates']
    rows[-1]['modeled_hourly_allocation_USD'] = rows[-1]['search_allocation_USD'] * 3600 / rows[-1]['search_wall_seconds']
    rows[-1]['storage_expansion_receipt_sha256'] = hashlib.sha256((nr / 'nr-volume-expansion-20261022.json').read_bytes()).hexdigest()
    lookup = {x['method']: x for x in rows}
    result = dict(searches=rows, observed_tblastx_to_DIAMOND_per_query_cost_ratio=lookup['tblastx']['search_allocation_USD_per_query'] / lookup['DIAMOND-sensitive']['search_allocation_USD_per_query'], observed_tblastx_to_blastn_total_cost_ratio=lookup['tblastx']['search_allocation_USD'] / lookup['blastn']['search_allocation_USD'], source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), interpretation=['Measured completed search wall time times allocated EC2, EBS, and IPv4 rates. The NR disk expansion is charged only for its overlap with the recorded search interval; the displayed effective hourly rate is duration weighted. Same conservative720hour month as the lifetime ledger.', 'Search allocation is part of lifetime spend, not an additional charge. Failed attempts, reference staging, idle time, final analysis, and backups are represented separately in the lifetime ledger.', 'NR and core_nt use different references, algorithms, query counts, thresholds, and hardware. The observed per-query ratio is not an isolated algorithm-speed comparison or a prediction for equal-sized batches.', 'The NR denominator includes all817744 pooled query mates, including other previously frozen experimental subsets and controls; no proportional allocation is presented as an independently timed200sample search.', 'Primary search memory includes mapped reference pages. Peak RSS does not establish a minimum RAM requirement.', 'No persistent translated core_nt index was built. Reference preparation/index/storage costs require their separate recorded inventories and timings; they are not included in search wall time.', 'The pilot uses a shared existing reference cache. A brief CPU profile and tiny characterization jobs overlapped its TBLASTX run.'])
    Path(a.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['searches', 'interpretation']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['pilot', 'nr', 'output']:
        p.add_argument('--' + name, required=True)
    main(p.parse_args())

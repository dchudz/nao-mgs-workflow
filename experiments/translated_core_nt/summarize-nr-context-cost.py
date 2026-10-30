"""Describe completed query-context reproducibility and allocate measured search cost."""
import datetime, hashlib, json, re
from pathlib import Path

root = Path('/home/ec2-user/experiment/protein')
small = root / 'seed-nr-standalone'
assert (small / 'complete.json').exists()
out = small / 'context-cost-and-mode.json'
assert not out.exists()
expansion = datetime.datetime.fromisoformat('2026-10-22T19:41:42+00:00')
records = {}
commands = []
for name in ('combined-validation', 'seed-nr-standalone'):
    d = root / name
    protocol = json.loads((d / 'search-protocol.json').read_text())
    log = (d / 'nr-search.log').read_text()
    command = protocol['command'].copy()
    for flag in ('--query', '--out'):
        command[command.index(flag)+1] = '<path>'
    commands.append(command)
    value = re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)', log)[-1]
    wall = sum(float(x)*60**i for i,x in enumerate(value.split(':')[::-1]))
    start = datetime.datetime.fromisoformat(protocol['started_utc'])
    end = start + datetime.timedelta(seconds=wall)
    before = max(0, min(wall, (expansion-start).total_seconds()))
    after = wall-before
    rates = {gib: 2.4344 + (gib*.08+3000*.005+875*.04)/720 + .005 for gib in (1280,1536)}
    allocation = (before*rates[1280] + after*rates[1536])/3600
    stages = json.loads((d / 'stage-timings.json').read_text())
    records[name] = dict(
        queries=protocol['queries'], query_sha256=protocol['query_sha256'],
        started_utc=protocol['started_utc'], GNU_wall_seconds=wall,
        GNU_user_seconds=float(re.findall(r'User time \(seconds\):\s*([\d.]+)',log)[-1]),
        GNU_system_seconds=float(re.findall(r'System time \(seconds\):\s*([\d.]+)',log)[-1]),
        GNU_max_RSS_KiB=int(re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)',log)[-1]),
        algorithms=sorted(set(re.findall(r'^Algorithm:\s*(.*)$',log,re.M))),
        allocated_compute_EBS_IPv4_USD=allocation,
        allocated_USD_per_query=allocation/protocol['queries'],
        search_seconds_before_volume_expansion=before, search_seconds_after_volume_expansion=after,
        logged_stage_groups=stages['groups'],
        input_hashes={n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in ['search-protocol.json','nr-search.log','stage-timings.json']},
    )
assert commands[0] == commands[1], 'Search options differ beyond input and output paths'
comparison = json.loads((small / 'context-comparison.json').read_text())
result = dict(
    generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    commands_identical_except_query_and_output_paths=True,
    executable_sha256=hashlib.sha256(Path(commands[0][0]).read_bytes()).hexdigest(),
    searches=records,
    reproducibility={k:v for k,v in comparison.items() if k!='raw_alignment_differing_query_IDs'},
    rates_USD_per_hour=rates,
    limits=[
        'The same executable options automatically selected different search algorithms as query context changed. This experiment does not isolate a batch-size effect from that automatic mode choice.',
        'Each cost is an allocation of GNU search elapsed time on r7a.8xlarge, plus provisioned gp3 and public IPv4. It excludes analysis, backup, staging, idle time, S3, requests, transfer and taxes.',
        'Use the experiment lifetime ledger for cap accounting; do not add these search allocations to that ledger.',
        'Original protocol rates used 1280GiB, a730-hour month and omitted IPv4. This supplementary record applies the recorded1280-to1536GiB expansion and a conservative720-hour month without editing those protocols.',
        'Completed_utc includes downstream analysis and upload; GNU elapsed search time is used instead.',
        'Logged stage durations are elapsed stage measurements, not CPU time, and may not sum exactly to GNU elapsed time.',
        'No persistent protein index was built for these searches. Reference seed structures were rebuilt during each search.',
        'Pooled primary classifications remain the reported primary evidence; standalone results are a separate reproducibility measurement.',
    ],
)
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'searches':{n:{k:v for k,v in r.items() if k not in ['input_hashes','logged_stage_groups']} for n,r in records.items()},'reproducibility':result['reproducibility']},indent=2))

"""One bounded mode control on the exact completed NR query subset; preserve every prior result."""
import argparse, collections, csv, datetime, hashlib, json, os, re, signal, subprocess, time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--commit', required=True)
a = p.parse_args()
assert re.fullmatch('[0-9a-f]{40}', a.commit)
root = Path('/home/ec2-user/experiment')
protein = root / 'protein'
auto = protein / 'seed-nr-standalone'
pooled = protein / 'combined-validation'
assert (auto / 'complete.json').exists()
assert not subprocess.run(['pgrep', '-x', 'diamond'], capture_output=True).stdout.strip(), 'Another DIAMOND process is active'
assert subprocess.check_output(['aws','--profile','pipeline-codex','--region','us-east-1','sts','get-caller-identity','--query','Account','--output','text'],text=True).strip() == '839993661825'
d = protein / 'seed-nr-forced-double-index'
d.mkdir(exist_ok=False)
for name in ('queries.fna', 'seed-selection.tsv'):
    (d / name).write_bytes((auto / name).read_bytes())
query_hash = hashlib.sha256((d / 'queries.fna').read_bytes()).hexdigest()
assert query_hash == 'b5d74d481f0291a671af931ff0e14f13b97f58f38ecd85b2bd2950e6cb3b841a'
cmd = json.loads((auto / 'search-protocol.json').read_text())['command'].copy()
assert '--algo' not in cmd
cmd[cmd.index('--query')+1] = str(d / 'queries.fna')
cmd[cmd.index('--out')+1] = str(d / 'nr.tsv')
cmd += ['--algo', '0']
start_utc = datetime.datetime.now(datetime.timezone.utc)
maximum = min(25200, (datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00')-start_utc).total_seconds())
assert maximum > 600
protocol = dict(
    source_commit=a.commit, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    command=cmd, started_utc=start_utc.isoformat(), queries=1436, pairs=746,
    query_sha256=query_hash, maximum_seconds=maximum, maximum_search_allocation_USD=maximum/3600*2.679511111111111,
    purpose='Separate the observed automatic query-indexed versus double-indexed mode difference from the surrounding-query-pool difference.',
    trigger='Completed same-options standalone repeat differs from pooled results in47of746pair statuses and74taxonomic assignments.',
    change='Only --algo0 is added to the completed standalone search options; query bytes, full NR snapshot, executable,32threads, masking, scores and output filters are retained.',
    limits=['One existing-query reproducibility control, not a larger NR biological panel or a replacement for same-core_nt work.', 'Prior pooled and standalone results remain unchanged.', 'If time-limited, preserve partial output and report the time bound; do not classify it as a completed search.'],
    option_source='https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options',
)
(d / 'search-protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
print(json.dumps(protocol), flush=True)
start = time.monotonic()
status = 'complete'
with (d / 'nr-search.log').open('w') as log:
    proc = subprocess.Popen(['/usr/bin/time', '-v']+cmd, stdout=log, stderr=log, start_new_session=True)
    try:
        code = proc.wait(timeout=maximum)
        if code: status = 'failed'
    except subprocess.TimeoutExpired:
        status = 'time_limit'
        os.killpg(proc.pid, signal.SIGTERM)
        try: proc.wait(timeout=15)
        except subprocess.TimeoutExpired: os.killpg(proc.pid, signal.SIGKILL); proc.wait()
        code = proc.returncode
    except BaseException:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
            try: proc.wait(timeout=15)
            except subprocess.TimeoutExpired: os.killpg(proc.pid, signal.SIGKILL); proc.wait()
        raise
resources = dict(status=status, returncode=code, wall_seconds=time.monotonic()-start, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(d / 'resources.json').write_text(json.dumps(resources, indent=2)+'\n')
print(json.dumps(resources), flush=True)
if status != 'complete': raise SystemExit(1)
assert 'Algorithm: Double-indexed' in (d / 'nr-search.log').read_text()
subprocess.run(['python3', str(root / 'summarize-seed-nr.py'), str(d)], check=True)
subprocess.run(['python3', str(root / 'summarize-diamond-stages.py'), str(d / 'nr-search.log'), str(d / 'stage-timings.json')], check=True)

def alignments(path):
    with path.open() as f: return collections.Counter(line.rstrip('\n') for line in f)

def calls(path):
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    result = {x['query_pair']:x for x in rows}
    assert len(result) == len(rows) == 746
    return result

new = alignments(d / 'nr.tsv')
after = calls(d / 'seed-validation.tsv')
comparisons = {}
changes = []
fields = ['nr_status', 'nr_taxid_lca', 'aligned_mates', 'taxonomic_relation']
for label, directory, raw in [('pooled_double', pooled, 'seed-alignments.tsv'), ('standalone_auto_query', auto, 'nr.tsv')]:
    old = alignments(directory / raw)
    before = calls(directory / 'seed-validation.tsv')
    assert before.keys() == after.keys()
    left_only = old-new
    right_only = new-old
    raw_queries = {line.split('\t',1)[0] for line in left_only.keys() | right_only.keys()}
    changed = [q for q in before if any(before[q][k] != after[q][k] for k in fields)]
    comparisons[label] = dict(
        alignment_multiset_equal=old==new, original_rows=sum(old.values()), fixed_mode_rows=sum(new.values()),
        rows_only_original=sum(left_only.values()), rows_only_fixed_mode=sum(right_only.values()),
        raw_differing_query_mates=len(raw_queries), pairs_with_summary_difference=len(changed),
        status_disagreements=sum(before[q]['nr_status']!=after[q]['nr_status'] for q in before),
        taxid_disagreements=sum(before[q]['nr_taxid_lca']!=after[q]['nr_taxid_lca'] for q in before),
        changed_status_transitions=dict(collections.Counter(before[q]['nr_status']+' -> '+after[q]['nr_status'] for q in before if before[q]['nr_status']!=after[q]['nr_status'])),
    )
    for q in changed:
        changes.append(dict(comparison=label, query_pair=q, **{'original_'+k:before[q][k] for k in fields}, **{'fixed_'+k:after[q][k] for k in fields}))
    del old, left_only, right_only
with (d / 'context-disagreements.tsv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=['comparison','query_pair']+['original_'+k for k in fields]+['fixed_'+k for k in fields], delimiter='\t'); w.writeheader(); w.writerows(changes)
(d / 'context-comparison.json').write_text(json.dumps(comparisons, indent=2)+'\n')
(d / 'complete.json').write_text(json.dumps(dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='complete', queries=1436, pairs=746), indent=2)+'\n')
print(json.dumps(comparisons, indent=2), flush=True)

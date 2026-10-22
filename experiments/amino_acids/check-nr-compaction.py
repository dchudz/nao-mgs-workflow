"""Check eligibility, scores, statuses and raw subset preservation against the original classifier."""
import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--root', default='/home/ec2-user/experiment')
p.add_argument('--binary', required=True)
p.add_argument('--output', required=True)
a = p.parse_args()
r, out = Path(a.root), Path(a.output)
out.mkdir(parents=True, exist_ok=False)
source = r / 'protein-nr-summary.py'
spec = importlib.util.spec_from_file_location('original', source)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
new_source=r/'protein-nr-summary.compact.py'
new_spec=importlib.util.spec_from_file_location('candidate',new_source)
n=importlib.util.module_from_spec(new_spec)
new_spec.loader.exec_module(n)

def rows(path):
    with open(path) as f:
        yield from csv.reader(f, delimiter='\t')

fixture = out / 'fixture.tsv'
fixture.write_text(
    'p/1\ta\t100\t1e-20\t25\t80\t11;12\n'
    'p/1\tb\t95\t1e-20\t26\t85\t11\n'
    'p/1\tc\t90\t1e-20\t25\t85\t13\n'
    'p/1\td\t89.9\t1e-20\t25\t85\t14\n'
    'p/1\te\t999\t1e-20\t24\t100\t21\n'
    'p/2\tf\t100\t1e-20\t25\t80\t21\n'
    'unknown/1\tg\t100\t1e-20\t25\t80\tN/A;;\n'
    'unknown/1\th\t90\t1e-20\t25\t80\t11;bad\n'
    'NR_CONTROL_c/1\ti\t50\t1e-20\t25\t80\t11;11\n'
    'SEED_s/2\tj\t55\t1e-20\t25\t80\t12\n'
    'short/1\tk\t1000\t1e-20\t24\t100\t11\n')
parents = {1: 1, 10239: 1, 11: 10239, 12: 10239, 13: 10239, 14: 10239, 21: 1}
with fixture.open('a') as f:
    for taxon in range(100000,125000): f.write(f'wide/1\ta{taxon}\t100\t1e-20\t25\t80\t{taxon}\n')
originals = {'wide': 11, 'p': 11, 'unknown': 11, 'NR_CONTROL_c': 11, 'SEED_s': 12, 'short': 11, 'absent': 11}

def check(path, dest, taxonomy=None, calls=None):
    subprocess.run([a.binary, str(path), str(dest), '0'], check=True)
    original = m.near_best_query_taxa(rows(path))
    with (dest / 'query-calls.tsv').open() as f:
        compact = {q: ({int(t) for t in taxa.split(';')},float(score)) for q,score,taxa in (line.rstrip('\n').split('\t') for line in f)}
    assert original == compact, 'Every read must retain exactly the same near-best taxon set and best score'
    if taxonomy is not None:
        assert m.summarize(rows(path), taxonomy, calls) == n.summarize_calls(compact, taxonomy, calls)
    for prefix, name in [('NR_CONTROL_', 'native-control-alignments.tsv'), ('SEED_', 'seed-alignments.tsv')]:
        h = hashlib.sha256()
        with path.open('rb') as f:
            for line in f:
                if line.startswith(prefix.encode()): h.update(line)
        assert h.hexdigest() == hashlib.sha256((dest / name).read_bytes()).hexdigest()
    return {'all_classification_records_identical': taxonomy is not None, 'classification_pairs_compared': len(calls) if calls is not None else 0, 'aligned_eligible_queries': len(original), 'counts': json.loads((dest / 'counts.json').read_text()), 'all_query_taxa_and_scores_identical': True, 'raw_control_and_seed_subset_bytes_identical': True}

synthetic = check(fixture, out / 'synthetic', parents, originals)
node_fixture=out/'nodes.dmp'
node_fixture.write_text(''.join(f'{t}\t|\t{parent}\t|\tno rank\t|\n' for t,parent in parents.items()))
call_fixture=out/'viral-calls.tsv'
call_fixture.write_text('query_pair\tprotein_taxid_lca\n'+''.join(f'{q}\t{t}\n' for q,t in originals.items()))
cli_output=out/'cli-classifications.tsv'
subprocess.run(['python3',str(new_source),'--nodes',str(node_fixture),'--query-calls',str(out/'synthetic/query-calls.tsv'),'--viral-calls',str(call_fixture),'--output',str(cli_output)],check=True,stdout=subprocess.DEVNULL)
with cli_output.open() as f: cli=list(csv.DictReader(f,delimiter='\t'))
assert cli==[{k:str(v) for k,v in row.items()} for row in m.summarize(rows(fixture),parents,originals)]
synthetic['CLI_large_taxon_field_classifications_identical']=True
bad = out / 'noncontiguous.tsv'
bad.write_text(fixture.read_text() + 'p/1\ta\t100\t1e-20\t25\t80\t11\n')
invalid = subprocess.run([a.binary, str(bad), str(out / 'rejected-noncontiguous'), '0'], capture_output=True, text=True)
assert invalid.returncode != 0 and 'contiguous' in invalid.stderr
real = r / 'protein/nr-thread-parity-memory-restart/threads32.tsv'
parents={}
for line in (r/'protein/current-nr/nodes.dmp').open():
    f=line.split('|');parents[int(f[0])]=int(f[1])
originals={}
for name in ['development-viral-pairs.tsv','full-viral-pairs.tsv']:
    with (r/'protein/combined-validation'/name).open() as f:
        for row in csv.DictReader(f,delimiter='\t'): originals[row['query_pair']]=int(row['protein_taxid_lca'])
actual = check(real, out / 'real-volume',parents,originals)
result = {'synthetic': synthetic, 'real_NR_volume': actual, 'noncontiguous_input_rejected': True, 'candidate_classifier_sha256': hashlib.sha256(new_source.read_bytes()).hexdigest(), 'original_classifier_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'binary_sha256': hashlib.sha256(Path(a.binary).read_bytes()).hexdigest(), 'limits': ['Query-call rows contain only the near-best taxon set and best eligible score, not synthetic alignments.', 'All raw alignments remain separate. Raw control and seed subsets preserve every input row.', 'Whole-input byte and expected alignment counts must pass before using a full compaction.']}
(out / 'parity.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))

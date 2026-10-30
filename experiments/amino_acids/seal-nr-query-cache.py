"""Seal a completed compaction only after full-input counts, parity and hashes pass."""
import argparse
import datetime
import json
import re
from pathlib import Path
from nr_evidence_inputs import digest

p = argparse.ArgumentParser()
p.add_argument('--directory', required=True)
p.add_argument('--cache', required=True)
p.add_argument('--parity', required=True)
p.add_argument('--commit', required=True)
a = p.parse_args()
assert re.fullmatch('[0-9a-f]{40}', a.commit)
d, cache = Path(a.directory), Path(a.cache)
assert not (d / 'query-call-cache.json').exists()
log = (d / 'nr-search.log').read_text()
assert re.search(r'Exit status:\s*0\s*$', log)
expected = int(re.search(r'Reported (\d+) pairwise alignments', log)[1])
queries = int(re.search(r'(\d+) queries aligned\.', log)[1])
counts = json.loads((cache / 'counts.json').read_text())
assert counts['input_rows'] == expected and counts['aligned_query_mates'] == queries
parity = json.loads(Path(a.parity).read_text())
assert parity['synthetic']['all_query_taxa_and_scores_identical']
assert parity['real_NR_volume']['all_query_taxa_and_scores_identical']
assert parity['synthetic']['all_classification_records_identical']
assert parity['real_NR_volume']['all_classification_records_identical']
assert parity['noncontiguous_input_rejected']
root = Path('/home/ec2-user/experiment')
assert digest(root / 'compact-nr-evidence') == parity['binary_sha256']
assert digest(root / 'protein-nr-summary.py') == parity['original_classifier_sha256']
assert digest(root / 'protein-nr-summary.compact.py') == parity['candidate_classifier_sha256']
raw = d / 'nr.tsv'
before = raw.stat()
assert before.st_size == counts['input_bytes']
raw_hash = digest(raw)
after = raw.stat()
assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
files = {}
for key, name in [('query_calls', 'query-calls.tsv'), ('native_controls', 'native-control-alignments.tsv'), ('seeds', 'seed-alignments.tsv')]:
    x = cache / name
    files[key] = dict(path=str(x), bytes=x.stat().st_size, sha256=digest(x))
result = dict(format='NR near-best query calls v1', completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_commit=a.commit, raw=dict(path=str(raw), bytes=before.st_size, mtime_ns=before.st_mtime_ns, sha256=raw_hash), counts=counts, files=files, parity=dict(path=str(Path(a.parity)), sha256=digest(a.parity), all_query_taxa_scores_and_classifications_identical=True, original_classifier_sha256=parity['original_classifier_sha256'], candidate_classifier_sha256=parity['candidate_classifier_sha256'], binary_sha256=parity['binary_sha256']), limits=['Raw search output is retained without alteration.', 'One query-call row contains the original eligible near-best taxon set and best score; it is not an alignment.', 'The original full-file Python aggregation continues independently for a whole-run comparison.'])
staged = d / 'query-call-cache.json.partial'
staged.write_text(json.dumps(result, indent=2) + '\n')
staged.replace(d / 'query-call-cache.json')
print(json.dumps(result))

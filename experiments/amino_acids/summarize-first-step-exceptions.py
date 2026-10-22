"""Summarize the completed raw/trimmed/reference audit without changing observed calls."""
import argparse, collections, json
from pathlib import Path
p = argparse.ArgumentParser(); p.add_argument('directory'); a = p.parse_args(); r = Path(a.directory)
boundary = json.loads((r / 'summary.json').read_text())
correction = json.loads((r / 'correction-diagnostic.json').read_text())['pairs']
references = json.loads((r / 'reference-masking-diagnostic.json').read_text())['pairs']
sequences = (r / 'trimmed-diagnostic.fastq').read_text().splitlines()[1::4]
assert all(set(sequence) <= set('ACGT') for sequence in sequences), 'Reference substring audit assumes unambiguous query bases'
by_id = collections.defaultdict(list)
for row in references:
    by_id[row['seq_id']].append(row)
assert {x['seq_id'] for x in correction} == set(by_id)
assert all(not x['raw_nt20_match'] for x in correction)
rest = [x for x in correction if not x['after_FASTP_nt20_match']]
changed = [x for x in correction if x['after_FASTP_nt20_match']]
assert all(max(x['minimum_substitutions_after_end_trimming_per_mate']) > 0 for x in changed)
masking = sum(any(max(x['alignment_DUST_masked_longest_exact_nt_per_mate']) >= 20 and max(x['screen_human_masked_longest_exact_nt_per_mate']) < 20 for x in by_id[row['seq_id']]) for row in rest)
assert boundary['added_beyond_nt20'] == len(correction) == boundary['stage']['screen']
summary = {'cohort_samples': boundary['completed_cohort_samples'], 'union_added_beyond_NT20': len(correction), 'samples': boundary['samples'], 'absent_from_raw_NT20': len(correction), 'gained_NT20_match_after_FASTP': len(changed), 'still_no_NT20_match_after_FASTP': len(rest), 'remaining_with_20nt_match_in_assigned_DUST_reference_but_not_human_masked_reference': masking, 'interpretation': 'This post hoc audit measures exact seed availability in existing raw and FASTP-processed queries against the frozen screening index and reported alignment references. End trimming cannot create a new substring match. Reference masking can change which exact seeds are available. These observations do not establish biological correctness or prospective sensitivity.'}
(r / 'mechanism-summary.json').write_text(json.dumps(summary, indent=2)); print(json.dumps(summary))

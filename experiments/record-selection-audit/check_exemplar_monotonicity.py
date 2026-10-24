"""Synthetic record-selection check; no biological reference or sample input."""
import collections, csv, gzip, hashlib, json, os, random, subprocess, time
from pathlib import Path

out = Path('/home/ec2-user/aux-exemplar-monotonicity-20261024-v2')
out.mkdir(exist_ok=False)
binary_root = Path('/home/ec2-user/experiment/rust-tools/target-musl/release')
rng = random.Random(284103)
fields = ['seq_id', 'prim_align_genome_id_all', 'prim_align_ref_start_unclipped',
          'prim_align_ref_start_unclipped_rev', 'prim_align_ref_end_unclipped',
          'prim_align_ref_end_unclipped_rev', 'prim_align_query_rc', 'prim_align_query_rc_rev',
          'query_seq', 'query_seq_rev', 'query_qual', 'query_qual_rev', 'synthetic_category']
rows = []
for group in range(3):
    sequences = [''.join(rng.choice('ACGT') for _ in range(160)) for _ in range(2)]
    for label, quality in [('unresolved', 'I'), ('class_X', '5')]:
        seqs = list(sequences)
        if label == 'class_X':
            seqs[0] = ('A' if seqs[0][0] != 'A' else 'C') + seqs[0][1:]
        rows.append(dict(zip(fields, [f'{group}_{label}', f'synthetic_group_{group}',
            '100', '300', '260', '460', 'False', 'True', *seqs,
            quality * 160, quality * 160, label])))

results = []
commands = []
for name, subset in [('all_records', rows), ('remove_unresolved_before_grouping', [r for r in rows if r['synthetic_category'] != 'unresolved'])]:
    source = out / (name + '.tsv.gz')
    with gzip.open(source, 'wt') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader(); writer.writerows(sorted(subset, key=lambda r: r['seq_id']))
    aligned, metadata = out / (name + '.aligned.tsv.gz'), out / (name + '.groups.tsv')
    final = out / (name + '.final.tsv.gz')
    cmds = [[str(binary_root / 'mark_duplicates'), '-i', str(source), '-o', str(aligned), '-m', str(metadata), '--num-threads', '1'],
            [str(binary_root / 'mark_duplicates_similarity'), '-i', str(aligned), '-o', str(final)]]
    start = time.monotonic()
    for i, command in enumerate(cmds):
        with (out / f'{name}.{i}.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=60,
                           env={**os.environ, 'RAYON_NUM_THREADS': '1'})
        commands.append(command)
    with gzip.open(final, 'rt') as f:
        selected = [r for r in csv.DictReader(f, delimiter='\t')
                    if r['seq_id'] == r['prim_align_dup_exemplar'] == r['sim_dup_exemplar']]
    counts = collections.Counter(r['synthetic_category'] for r in selected)
    results.append({'case': name, 'input_records': len(subset), 'joint_exemplars': len(selected),
                    'selected_categories': dict(counts), 'class_X_meets_three_exemplar_threshold': counts['class_X'] >= 3,
                    'wall_seconds': time.monotonic() - start})
record = {'data': 'Six synthetic records, random strings, invented categories and reference keys; no biological sample or reference used.',
          'results': results, 'commands': commands,
          'binary_sha256': {name: hashlib.sha256((binary_root / name).read_bytes()).hexdigest()
                            for name in ['mark_duplicates', 'mark_duplicates_similarity']},
          'scope': 'A category-support counterexample for duplicate exemplar selection. Does not establish that this transition occurs in the experimental biological sample.',
          'implication': 'Filtering before exemplar selection can change which category gets support. A pre-filter category threshold is not automatically an upper bound on the post-filter category threshold.'}
(out / 'result.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))

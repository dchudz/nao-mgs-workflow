"""Count species support after a common similarity-only duplicate check on both arms.

This deliberately does not fabricate protein alignment-duplicate coordinates or
claim to reproduce the production two-stage exemplar definition.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

IMAGE = 'public.ecr.aws/q0n1c7g8/nao-mgs-workflow/rust-tools@sha256:e69b93d9a1a0ae41c59ec4c0a6a0841f3a7f4871fc2f9b333eaa094563612270'


def main(a):
    r, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    helper = r / 'summarize-nr-species-calls.py'
    spec = importlib.util.spec_from_file_location('species_counts', helper)
    counts = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(counts)
    table = counts.table
    sys.path.insert(0, str(r / 'experiments/amino_acids'))
    from union_pairs import pairs
    samples = [x['run_accession'] for x in json.loads((r / 'heldout-ranges.json').read_text())][:200]
    assert len(samples) == len(set(samples)) == 200
    resolution = {x['query_pair']: x for x in table(a.resolution)}
    meta = {x['query_pair']: x for x in table(r / 'cohort-storage/protein200/all.tsv.gz') if x['query_pair'] in resolution}
    assert meta.keys() == resolution.keys()
    sample_queries = collections.defaultdict(list)
    for q, x in resolution.items():
        assert x['sample'] in samples
        sample_queries[x['sample']].append(q)
    members = collections.defaultdict(set)
    for x in table(r / 'cohort-storage/protein200/all.members.tsv.gz'):
        if x['query_pair'] in resolution:
            members[x['query_pair']].add(x['seq_id'])
    assert all(len(members[q]) == int(meta[q]['occurrences']) for q in meta)
    batches = {s: p.parent for p in (r / 'cohort-storage/batches').glob('*/complete.json') for s in json.loads(p.read_text())['samples']}
    parents, ranks = {}, {}
    for line in (r / 'reference/results/taxonomy-nodes.dmp').open():
        f = [x.strip() for x in line.split('|')]
        parents[int(f[0])], ranks[int(f[0])] = int(f[1]), f[2]
    def species(t):
        seen = set()
        while t in parents and t not in seen:
            if ranks[t] == 'species':
                return t
            seen.add(t)
            t = parents[t]
        return None
    human = {int(x['taxid']): x['infection_status_human'] for x in table(r / 'reference/results/total-virus-db-annotated.tsv.gz') if x['rank'] == 'species'}
    baseline = {}
    tiers = {name: {} for name in ['NR_resolved_species', 'NR_and_viral_reference_same_species']}
    metrics = []
    mapping = gzip.open(out / 'cluster-membership.tsv.gz', 'wt')
    mw = csv.writer(mapping, delimiter='\t')
    mw.writerow(['sample', 'canonical_pair', 'similarity_exemplar', 'native_species', 'protein_query', 'nr_species'])
    for sample in samples:
        start = time.monotonic()
        batch = batches[sample]
        paths = list((batch / 'output/marked').glob(f'*_{sample}_{a.native_arm}_duplicate_reads.tsv.gz'))
        assert len(paths) == 1
        native = {x['seq_id']: x for x in table(paths[0])}
        quality, native_species, seen_native = {}, collections.defaultdict(set), set()
        def retain(key):
            name, s1, q1, s2, q2 = [x.decode() for x in key]
            seqs = s1, s2
            candidate = sum(map(ord, q1 + q2)), q1, q2, name
            if seqs not in quality or candidate > quality[seqs]:
                quality[seqs] = candidate
            return seqs
        for key, _ in pairs(batch / 'output/fastp' / f'{sample}_{a.native_arm}_fastp.fastq.gz'):
            name = key[0].decode().lstrip('@')
            if name in native:
                seen_native.add(name)
                seqs = retain(key)
                t = species(int(native[name]['aligner_taxid_lca']))
                if t is not None:
                    native_species[seqs].add(t)
        assert seen_native == native.keys()
        query_seqs = {q: (meta[q]['query_seq'], meta[q]['query_seq_rev']) for q in sample_queries[sample]}
        wanted_ids = {name: q for q in query_seqs for name in members[q]}
        assert len(wanted_ids) == sum(len(members[q]) for q in query_seqs)
        seen_protein = set()
        for key, _ in pairs(batch / 'output/fastp' / f'{sample}_aa12_fastp.fastq.gz'):
            name = key[0].decode().lstrip('@')
            if name in wanted_ids:
                seqs = retain(key)
                assert seqs == query_seqs[wanted_ids[name]], (sample, name, 'Prepared sequences differ from source FASTP')
                seen_protein.add(name)
        assert seen_protein == wanted_ids.keys(), (sample, 'Missing original prepared-read qualities')
        ids = {seqs: 'p' + hashlib.sha256((seqs[0] + '\n' + seqs[1]).encode()).hexdigest() for seqs in quality}
        assert len(set(ids.values())) == len(ids)
        d = out / sample
        d.mkdir()
        with gzip.open(d / 'input.tsv.gz', 'wt') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['seq_id', 'query_seq', 'query_seq_rev', 'query_qual', 'query_qual_rev', 'prim_align_dup_exemplar'])
            for seqs in sorted(quality, key=ids.get):
                _, q1, q2, _ = quality[seqs]
                w.writerow([ids[seqs], *seqs, q1, q2, ids[seqs]])
        cmd = ['docker', 'run', '--rm', '--network', 'none', '--cpus', '1', '--memory', '8g', '-v', f'{r}:{r}', IMAGE, 'mark_duplicates_similarity', '-i', str(d / 'input.tsv.gz'), '-o', str(d / 'output.tsv.gz')]
        with (d / 'stderr.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=log, check=True, timeout=600)
        cluster = {x['seq_id']: x['sim_dup_exemplar'] for x in table(d / 'output.tsv.gz')}
        assert cluster.keys() == set(ids.values()) and set(cluster.values()) <= cluster.keys()
        b = collections.defaultdict(set)
        for seqs, taxa in native_species.items():
            for t in taxa:
                b[t].add(cluster[ids[seqs]])
        baseline[sample] = b
        for name in tiers:
            tiers[name][sample] = collections.defaultdict(set)
        protein_by_seq = {seqs: q for q, seqs in query_seqs.items()}
        assert len(protein_by_seq) == len(query_seqs)
        for q, seqs in query_seqs.items():
            x = resolution[q]
            if x['nr_species']:
                t = int(x['nr_species'])
                tiers['NR_resolved_species'][sample][t].add(cluster[ids[seqs]])
                if x['species_relation'] == 'same_species':
                    tiers['NR_and_viral_reference_same_species'][sample][t].add(cluster[ids[seqs]])
        for seqs, ident in ids.items():
            q = protein_by_seq.get(seqs, '')
            mw.writerow([sample, ident, cluster[ident], ';'.join(map(str, sorted(native_species[seqs]))), q, resolution[q]['nr_species'] if q else ''])
        metrics.append(dict(sample=sample, exact_pair_inputs=len(ids), similarity_clusters=len(set(cluster.values())), protein_viral_only_pairs=len(query_seqs), seconds=time.monotonic() - start))
        (d / 'complete.json').write_text(json.dumps(dict(metrics[-1], command=cmd), indent=2) + '\n')
        print(json.dumps(metrics[-1]), flush=True)
    mapping.close()
    rows, summary = counts.summarize(baseline, tiers, human)
    rows = [{k.replace('exact_pairs', 'similarity_clusters'): v for k, v in x.items()} for x in rows]
    with (out / 'sample-species-calls.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ['sample'], delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    result = dict(samples=len(samples), native_arm=a.native_arm, threshold_similarity_clusters=3, counts=summary, sample_metrics=metrics, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), helper_sha256=hashlib.sha256(helper.read_bytes()).hexdigest(), resolution_sha256=hashlib.sha256(Path(a.resolution).read_bytes()).hexdigest(), container=IMAGE, limits=['Supplementary common similarity-only support check; not the production alignment-plus-similarity exemplar definition.', 'Joint per-sample pool contains baseline retained sequences and every NR-viral-only candidate, including unresolved species. Each exact pair enters once, in deterministic sequence-hash order.', 'Original FASTP sequence and membership assertions recover real qualities; retain the maximum summed-quality observation per exact sequence pair, with lexical tie breaking.', 'Each read keeps its own taxonomic evidence. Cluster membership only reduces support counts and never transfers a validation verdict to another read.', 'A cluster can contain differently assigned reads. Per-species support counts do not establish a unique species assignment for the entire cluster.', 'Sequence support does not establish infection or known-truth accuracy.'])
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'sample_metrics'}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='/home/ec2-user/experiment')
    p.add_argument('--native-arm', choices=['nt', 'nt20'], default='nt20')
    p.add_argument('--resolution', required=True)
    p.add_argument('--output', required=True)
    main(p.parse_args())

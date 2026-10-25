"""Read frozen results; add species/exemplar tabulations without changing calls."""
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''): h.update(block)
    return h.hexdigest()


def read_table(path):
    with gzip.open(path, 'rt') as f:
        yield from csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)


def joint(row):
    return row['seq_id'] == row['prim_align_dup_exemplar'] == row['sim_dup_exemplar']


def transition(a, b, threshold):
    if a >= threshold: return 'retained' if b >= threshold else 'lost'
    if b >= threshold: return 'gained_from_zero' if a == 0 else 'promoted'
    return 'below_threshold_both'


def count_rows(rows, taxonomy):
    counts = collections.Counter(); unresolved = collections.Counter()
    total = exemplars = sim_only = 0
    for row in rows:
        total += 1
        sim_only += row['seq_id'] == row['sim_dup_exemplar']
        if not joint(row): continue
        exemplars += 1
        tax = taxonomy.get(row['aligner_taxid_lca'])
        assert tax is not None, row['aligner_taxid_lca']
        species = tax['taxid_species']
        if species not in ('', 'NA', '0', 'None'):
            assert taxonomy[species]['rank'] == 'species'
            counts[species] += 1
        else: unresolved[tax['rank']] += 1
    assert sum(counts.values()) + sum(unresolved.values()) == exemplars
    return dict(rows=total, joint_exemplars=exemplars, similarity_exemplars=sim_only,
                species_exemplars=dict(counts), unresolved_exemplars_by_rank=dict(unresolved))


def fastp_counts(d, arm):
    if arm == 'production':
        with (d/'baseline-run/trace.tsv').open() as f:
            tasks = [r for r in csv.DictReader(f, delimiter='\t')
                     if r['process'].endswith('EXTRACT_VIRAL_READS_SHORT:FASTP')]
        assert len(tasks) == 1
        paths = list(Path(tasks[0]['workdir']).glob('*.json'))
    else:
        paths = list((d/'prepare/published/FASTP').glob('*.json'))
    assert len(paths) == 1, paths
    data = json.loads(paths[0].read_text())['summary']
    before = data['before_filtering']['total_reads']; after = data['after_filtering']['total_reads']
    assert before % 2 == after % 2 == 0
    return dict(before_pairs=before//2, after_pairs=after//2, file=str(paths[0]), sha256=sha(paths[0]))


def main():
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    taxonomy={r['taxid']:r for r in read_table(a.reference)}; ref_hash=sha(a.reference)
    completed=[]
    for proof in sorted((a.root/'comparisons').glob('*/completed-evidence-backup.json')):
        d=proof.parent; target=a.out/d.name
        assert json.loads(proof.read_text())['readback_verified']
        if (target/'complete.json').exists(): completed.append(d.name); continue
        expected=json.loads((d/'comparison-report/validation-status-counts.json').read_text())
        tables={}; arms={}; filters={}
        for arm in ('production','protein_variant'):
            paths=list((d/(arm+'-downstream/output/results_downstream')).glob('*_validation_hits.tsv.gz'))
            assert len(paths)==1; tables[arm]=paths[0]
            arms[arm]=count_rows(read_table(paths[0]),taxonomy)
            assert all(arms[arm][k]==expected[arm][k] for k in ('rows','joint_exemplars'))
            filters[arm]=fastp_counts(d,arm)
        target.mkdir(exist_ok=True)
        species=set(arms['production']['species_exemplars'])|set(arms['protein_variant']['species_exemplars'])
        records=[]; gained=set()
        for species_id in sorted(species,key=int):
            x=arms['production']['species_exemplars'].get(species_id,0)
            y=arms['protein_variant']['species_exemplars'].get(species_id,0)
            rec=dict(group=d.name,species_taxid=species_id,species_name=taxonomy[species_id]['name'],
                     production_exemplars=x,protein_variant_exemplars=y)
            for n in (2,3,5): rec['threshold_'+str(n)]=transition(x,y,n)
            records.append(rec)
            if rec['threshold_3']=='gained_from_zero': gained.add(species_id)
        fields=['group','species_taxid','species_name','production_exemplars','protein_variant_exemplars','threshold_2','threshold_3','threshold_5']
        with (target/'species.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader();w.writerows(records)
        # These are the original, post-trim DNA observations for a fixed NR audit.
        query_count=pair_count=0
        with gzip.open(target/'gained-exemplars.tsv.gz','wt') as out:
            fields=['group','seq_id','species_taxid','query_seq','query_seq_rev']
            w=csv.DictWriter(out,fieldnames=fields,delimiter='\t');w.writeheader()
            for row in read_table(tables['protein_variant']):
                if not joint(row):continue
                sp=taxonomy[row['aligner_taxid_lca']]['taxid_species']
                if sp not in gained:continue
                w.writerow(dict(group=d.name,seq_id=row['seq_id'],species_taxid=sp,
                                query_seq=row['query_seq'],query_seq_rev=row['query_seq_rev']))
                pair_count+=1
                query_count+=sum(row[k] not in ('','NA') for k in ('query_seq','query_seq_rev'))
        assert pair_count==sum(arms['protein_variant']['species_exemplars'][s] for s in gained)
        receipt=dict(label=d.name,created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                     pipeline_source_commit=json.loads((d/'complete.json').read_text())['source_commit'],
                     analysis_sha256=sha(Path(__file__)),reference_sha256=ref_hash,
                     input_tables={arm:dict(path=str(path),sha256=sha(path)) for arm,path in tables.items()},
                     arms=arms,fastp=filters,gained_from_zero_threshold_3_calls=len(gained),
                     gained_exemplar_pairs=pair_count,gained_mate_queries=query_count,
                     thresholds=[2,3,5],species_rule='Map aligner_taxid_lca to reference taxid_species; never force an above-species assignment to species.',
                     calling_rule='Additional reporting threshold, not a binary output emitted by the frozen workflow.',
                     outputs={name:sha(target/name) for name in ('species.tsv','gained-exemplars.tsv.gz')})
        (target/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n');completed.append(d.name)
        print(json.dumps({'label':d.name,'species':len(species),'gained_calls':len(gained),'exemplar_pairs':pair_count,'mate_queries':query_count}),flush=True)
    print(json.dumps({'completed_labels':completed}),flush=True)


if __name__=='__main__':main()

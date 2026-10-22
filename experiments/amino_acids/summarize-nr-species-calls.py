"""Supplement original protein summaries with species calls supported by NR resolution itself."""
import argparse,collections,csv,gzip,hashlib,json,sys
from pathlib import Path

def table(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:yield from csv.DictReader(f,delimiter='\t')
def summarize(baseline,tiers,human):
    records=[];counts={}
    for tier,support in tiers.items():
        c=collections.Counter()
        for sample,base in baseline.items():
            called=support.get(sample,{})
            for species in sorted(base.keys()|called.keys()):
                b=base.get(species,set());p=called.get(species,set());before=len(b)>=3;after=len(p)>=3;union=len(b|p)>=3
                if before or after or union:
                    row=dict(sample=sample,tier=tier,species=species,human_annotation=human.get(species,'NA'),baseline_exact_pairs=len(b),protein_exact_pairs=len(p),union_exact_pairs=len(b|p),baseline_call=before,protein_call=after,union_call=union,gained_protein_call=after and not before,gained_union_call=union and not before,lost_protein_call=before and not after)
                    records.append(row)
                    for key in ['baseline_call','protein_call','union_call','gained_protein_call','gained_union_call','lost_protein_call']:
                        c[key]+=int(row[key]);c['human_annotated_'+key]+=int(row[key] and human.get(species)=='1')
                    if row['gained_protein_call']:c['gained_protein_call_baseline_support_'+str(len(b))]+=1
                    if row['gained_union_call']:c['gained_union_call_baseline_support_'+str(len(b))]+=1
        counts[tier]=dict(c)
    return records,counts

def main(a):
    r=Path(a.root);d=r/'cohort-storage/protein200';out=Path(a.output);out.mkdir(parents=True,exist_ok=False);sys.path.insert(0,str(r/'experiments/amino_acids'));from union_pairs import pairs
    evidence=list(table(a.evidence));resolution={x['query_pair']:x for x in table(a.resolution)};viral={x['query_pair']:x for x in evidence if x['nr_status']=='viral_only'};assert viral.keys()==resolution.keys(),'NR-resolution rows must cover every viral-only query pair exactly'
    samples=[x['run_accession'] for x in json.loads((r/'heldout-ranges.json').read_text())][:200];assert len(samples)==len(set(samples))==200
    meta={x['query_pair']:(x['query_seq'],x['query_seq_rev']) for x in table(d/'all.tsv.gz') if x['query_pair'] in viral};assert meta.keys()==viral.keys()
    batches={s:m.parent for m in (r/'cohort-storage/batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']};assert set(samples)<=batches.keys()
    parents={};ranks={}
    for line in (r/'reference/results/taxonomy-nodes.dmp').open():
        f=[x.strip() for x in line.split('|')];parents[int(f[0])]=int(f[1]);ranks[int(f[0])]=f[2]
    def species(t):
        seen=set()
        while t in parents and t not in seen:
            if ranks[t]=='species':return t
            seen.add(t);t=parents[t]
        return None
    annotations={int(x['taxid']):x for x in table(r/'reference/results/total-virus-db-annotated.tsv.gz') if x['rank']=='species'};human={t:x['infection_status_human'] for t,x in annotations.items()}
    tiers={name:collections.defaultdict(lambda:collections.defaultdict(set)) for name in ['NR_resolved_species','NR_and_viral_reference_same_species']}
    for q,x in resolution.items():
        if not x['nr_species']:continue
        t=int(x['nr_species']);sample=x['sample'];assert sample in samples and sample==viral[q]['sample'];tiers['NR_resolved_species'][sample][t].add(meta[q])
        if x['species_relation']=='same_species':tiers['NR_and_viral_reference_same_species'][sample][t].add(meta[q])
    reports={};allrows=[]
    for arm in ['nt','nt20']:
        baseline={}
        for sample in samples:
            batch=batches[sample];paths=list((batch/'output/marked').glob(f'*_{sample}_{arm}_duplicate_reads.tsv.gz'));assert len(paths)==1
            native={x['seq_id']:x for x in table(paths[0])};calls=collections.defaultdict(set);seen=set()
            for key,raw in pairs(batch/'output/fastp'/f'{sample}_{arm}_fastp.fastq.gz'):
                q=key[0].decode().lstrip('@')
                if q not in native:continue
                seen.add(q);t=species(int(native[q]['aligner_taxid_lca']))
                if t is not None:calls[t].add((key[1].decode(),key[3].decode()))
            assert seen==native.keys();baseline[sample]=calls
        rows,counts=summarize(baseline,tiers,human);allrows.extend(dict(baseline_arm=arm,**x) for x in rows);reports[arm]=counts
    with (out/'sample-species-calls.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]) if allrows else ['baseline_arm','sample','tier','species'],delimiter='\t');w.writeheader();w.writerows(allrows)
    result=dict(samples=200,threshold_exact_pairs=3,viral_only_query_pairs=len(viral),viral_only_with_NR_species=sum(bool(x['nr_species']) for x in resolution.values()),viral_only_same_species=sum(x['species_relation']=='same_species' for x in resolution.values()),counts=reports,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),evidence_sha256=hashlib.sha256(Path(a.evidence).read_bytes()).hexdigest(),resolution_sha256=hashlib.sha256(Path(a.resolution).read_bytes()).hexdigest(),limits=['Supplementary analysis; original viral-reference-species summaries remain unchanged','NR-above-species LCAs do not count as species-resolved support','Common exact complete trimmed nucleotide-pair keys on both sides; these are not the production alignment/similarity duplicate exemplars','NR species mapping uses the frozen NR resolution table; native species mapping uses the original nucleotide-reference taxonomy','Protein calls and union-with-baseline calls are separate, explicitly labeled constructions','Sequence support and taxonomic agreement do not establish infection or known-truth accuracy'])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='/home/ec2-user/experiment');p.add_argument('--evidence',required=True);p.add_argument('--resolution',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(a)

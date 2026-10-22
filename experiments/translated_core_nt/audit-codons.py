"""Reconstruct translated hit sequences from the original nucleotide reference."""
import argparse,csv,importlib.util,json,subprocess
from pathlib import Path
from Bio.Seq import Seq
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def audit(row,queries,database):
    q=queries[row['qseqid']];a,b=int(row['qstart']),int(row['qend']);qnt=q[min(a,b)-1:max(a,b)]
    if int(row['qframe'])<0:qnt=str(Seq(qnt).reverse_complement())
    start,end=int(row['sstart']),int(row['send']);strand='minus' if int(row['sframe'])<0 else 'plus'
    command=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}:ro',m.IMAGE,'blastdbcmd','-db',str(database),'-entry',row['sseqid'],'-range',str(min(start,end))+'-'+str(max(start,end)),'-strand',strand,'-outfmt','%s']
    snt=''.join(subprocess.check_output(command,text=True).splitlines()).upper()
    qaa=str(Seq(qnt).translate(table=1));saa=str(Seq(snt).translate(table=1));assert qaa.upper()==row['qseq'].upper() and saa.upper()==row['sseq'].upper(),'Reference/query reconstruction disagrees with alignment'
    assert len(qnt)==len(snt)==3*int(row['length']) and '-' not in row['qseq']+row['sseq']
    same_aa_different_codon=sum(qaa[i]==saa[i] and qaa[i] in 'ACDEFGHIKLMNPQRSTVWY' and qnt[3*i:3*i+3]!=snt[3*i:3*i+3] for i in range(len(qaa)))
    longest=run=0
    for x,y in zip(qnt,snt):
        run=run+1 if x==y and x in 'ACGT' else 0;longest=max(longest,run)
    return dict(query=row['qseqid'],subject=row['sseqid'],query_frame=int(row['qframe']),subject_frame=int(row['sframe']),taxids=row['staxids'],bitscore=float(row['bitscore']),evalue=float(row['evalue']),alignment_aa=int(row['length']),reported_aa_identity_pct=float(row['pident']),aligned_nt_identity_pct=100*sum(x==y for x,y in zip(qnt,snt))/len(qnt),synonymous_changed_codons=same_aa_different_codon,longest_exact_nt_run_in_alignment=longest,coordinate_query_coverage_pct=100*len(qnt)/len(q),query_nt=qnt,subject_nt=snt,query_aa=qaa,subject_aa=saa,extraction_command=command)

def characterize():
    f=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=f/'codon-families';queries=m.fasta(d/'queries.fna')
    rows=[dict(zip(m.FIELDS,x)) for x in csv.reader((d/'tblastx.tsv').open(),delimiter='\t')]
    database=d/'audit-db'
    subprocess.run(['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE,'makeblastdb','-in',str(f/'reference.fna'),'-dbtype','nucl','-parse_seqids','-out',str(database)],check=True,stdout=subprocess.DEVNULL)
    all_results=[audit(x,queries,database) for x in rows]
    assert any(x['subject_frame']<0 for x in all_results)
    results=[x for x in all_results if x['query']!='reference_control' and x['alignment_aa']==80]
    assert len(results)==2 and all(x['synonymous_changed_codons']>0 and x['aligned_nt_identity_pct']==60 for x in results)
    assert {x['query_frame'] for x in results}=={1,-3}
    m.save(d/'codon-audit-characterization.json',dict(created_utc=m.now(),rows=all_results,reference_sha256=m.sha(f/'reference.fna'),excluded_from_real_data=True));print('CODON_RECONSTRUCTION_VERIFIED')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--characterize',action='store_true');p.add_argument('--directory',type=Path);p.add_argument('--analysis',default='analysis');p.add_argument('--method',choices=['tblastx','query-tblastn-fast'],default='tblastx');p.add_argument('--baseline',choices=['blastn','tblastx'],default='blastn');a=p.parse_args()
    if a.characterize:characterize()
    else:
        d=a.directory;assert d is not None;out=d/a.analysis;calls=list(csv.DictReader((out/'per-mate.tsv').open(),delimiter='\t'));lookup={(x['query'],x['method'],x['threshold']):x for x in calls}
        queries=m.fasta(d/'queries.fna');selected=[];changed=[];assert a.method!=a.baseline
        for q in queries:
            left=lookup[q,a.baseline,'strict_E1e-5'];right=lookup[q,a.method,'strict_E1e-5']
            if (left['status'],left['pipeline_natural_lca'])!=(right['status'],right['pipeline_natural_lca']):
                changed.append(q)
                if right['status']!='no_alignment':selected.append(q)
        rows=[x for x in csv.DictReader((out/'near-best-alignments.tsv').open(),delimiter='\t') if x['method']==a.method and x['threshold']=='strict_E1e-5' and x['qseqid'] in selected]
        # One highest-scoring subject per changed mate; all changed mate identities retained.
        chosen=[min((x for x in rows if x['qseqid']==q),key=lambda x:(-float(x['bitscore']),float(x['evalue']),x['sseqid'])) for q in sorted(selected)]
        results=[audit(x,queries,m.DB) for x in chosen]
        name='codon-audit.json' if (a.method,a.baseline)==('tblastx','blastn') else 'codon-audit-'+a.method+'-vs-'+a.baseline+'.json'
        m.save(out/name,dict(created_utc=m.now(),audited_method=a.method,baseline_method=a.baseline,all_disagreeing_mates=changed,changed_without_audited_method_hit=sorted(set(changed)-set(selected)),selection='Every mate with a strict audited-method hit whose status or natural LCA differs from the baseline; one deterministic highest-score subject per mate',rows=results,limits=['Alignment reconstruction explains observed sequence differences; it does not establish infection or a correct species assignment','No exact28nt run within an alignment alone is not proof of the sole reason a method missed the read']))
        print(json.dumps({'audited_mates':len(results),'with_synonymous_codon_changes':sum(x['synonymous_changed_codons']>0 for x in results)}))

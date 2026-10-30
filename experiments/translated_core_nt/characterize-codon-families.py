"""Test synonymous changes beyond third-position wobble and reverse-frame coordinates."""
import csv,importlib.util,json,subprocess
from pathlib import Path
from Bio.Seq import Seq
from Bio.Data import CodonTable
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
fixture=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=fixture/'codon-families';d.mkdir(exist_ok=False)
reference=next(iter(m.fasta(fixture/'reference.fna').values()));table=CodonTable.unambiguous_dna_by_id[1].forward_table;changed=[]
for i in range(0,len(reference),3):
    codon=reference[i:i+3];aa=table[codon];alternatives=[c for c,a in table.items() if a==aa]
    changed.append(max(alternatives,key=lambda c:(sum(a!=b for a,b in zip(c,codon)),c)))
recoded=''.join(changed);assert str(Seq(reference).translate())==str(Seq(recoded).translate())
queries={'reference_control':reference,'maximal_synonymous':recoded,'reverse_frame3_synonymous':str(Seq('AC'+recoded).reverse_complement())}
(d/'queries.fna').write_text(''.join('>'+q+'\n'+s+'\n' for q,s in queries.items()))
identity=100*sum(a==b for a,b in zip(reference,recoded))/len(reference);results={}
for name in ['megablast','blastn','dc-megablast','tblastx']:
    config=json.loads(Path(__file__).with_name(name+'.json').read_text())
    cmd=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE,config['executable']]+config['arguments']+['-db',str(fixture/'db'),'-query',str(d/'queries.fna'),'-outfmt',m.OUTFMT,'-qcov_hsp_perc','30','-evalue','10','-max_hsps','5','-num_alignments','250']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=60);assert result.returncode==0
    (d/(name+'.tsv')).write_text(result.stdout);rows=[dict(zip(m.FIELDS,x.split('\t'))) for x in result.stdout.splitlines()]
    groups={q:dict(rows=sum(x['qseqid']==q for x in rows),best_evalue=min((float(x['evalue']) for x in rows if x['qseqid']==q),default=None),max_query_HSP_coverage=max((float(x['qcovhsp']) for x in rows if x['qseqid']==q),default=None)) for q in queries}
    if name=='tblastx':
        for q in ['maximal_synonymous','reverse_frame3_synonymous']:
            exact=[x for x in rows if x['qseqid']==q and int(x['length'])==80 and float(x['pident'])==100]
            assert exact
            if q.startswith('reverse'):assert any(x['qframe']=='-3' and x['qstart']=='240' and x['qend']=='1' for x in exact)
    results[name]=dict(command=cmd,stderr=result.stderr,queries=groups)
m.save(d/'characterization.json',dict(created_utc=m.now(),reference_nt=len(reference),synonymous_nt_identity_pct=identity,protein_identity_pct=100,query_sha256=m.sha(d/'queries.fna'),results=results,interpretation='Tiny reference tests mechanisms and coordinate handling only; Evalues are not evidence about significance against comprehensive core_nt; no runtime claim'))
print(json.dumps(dict(synonymous_nt_identity_pct=identity,protein_identity_pct=100,methods={k:v['queries'] for k,v in results.items()})))

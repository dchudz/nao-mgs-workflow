"""Controlled synthetic sensitivity ladder; no population accuracy or full-reference significance claim."""
import csv,importlib.util,json,random,subprocess
from pathlib import Path
from Bio.Seq import Seq
from Bio.Data import CodonTable
spec=importlib.util.spec_from_file_location('translation',Path(__file__).with_name('characterize-query-translation.py'));t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);m=t.m
f=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=f/'divergence';d.mkdir(exist_ok=False)
reference=next(iter(m.fasta(f/'reference.fna').values()));protein=str(Seq(reference).translate());table=CodonTable.unambiguous_dna_by_id[1].forward_table;alphabet=sorted(set(table.values()));records={};truth={}
def encode(peptide):
    return ''.join(max((c for c,a in table.items() if a==aa),key=lambda c:(sum(x!=y for x,y in zip(c,reference[3*i:3*i+3])),c)) for i,aa in enumerate(peptide))
for count in [0,4,8,16,24,32]:
    for placement in ['even','random']:
        rng=random.Random(741203+count);positions=set((int((i+.5)*len(protein)/count) for i in range(count))) if count and placement=='even' else set(rng.sample(range(len(protein)),count))
        peptide=''.join(rng.choice([x for x in alphabet if x!=aa]) if i in positions else aa for i,aa in enumerate(protein));dna=encode(peptide);q=placement+'_'+str(count)+'substitutions';records[q]=dna;truth[q]=dict(AA_substitutions=count,AA_identity_pct=100*(1-count/len(protein)),nt_identity_pct=100*sum(a==b for a,b in zip(dna,reference))/len(reference),placement=placement)
rng=random.Random(741203);shuffled=list(protein);rng.shuffle(shuffled);records['shuffled_control']=encode(''.join(shuffled));truth['shuffled_control']=dict(control='Composition-preserving protein shuffle; not a biological negative truth label')
(d/'queries.fna').write_text(''.join('>'+q+'\n'+s+'\n' for q,s in records.items()));_,mapping=t.translate_queries(d/'queries.fna',d/'queries.faa');m.save(d/'query-map.json',mapping)
prefix=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE];results={}
methods=[json.loads(Path(__file__).with_name(name+'.json').read_text()) for name in ['megablast','blastn','dc-megablast','tblastx','query-tblastn-fast']]
methods += [dict(name='tblastx-word4',executable='tblastx',arguments=['-word_size','4','-threshold','18','-query_gencode','1','-db_gencode','1']),dict(methods[-1],name='query-tblastn-standard',arguments=['-task','tblastn']+methods[-1]['arguments'][2:])]
for method in methods:
    translated_query=method['name'].startswith('query-');query=d/('queries.faa' if translated_query else 'queries.fna')
    command=prefix+[method['executable']]+method['arguments']+['-db',str(f/'db'),'-query',str(query),'-outfmt',m.OUTFMT,'-qcov_hsp_perc','30','-evalue','10','-max_hsps','5','-num_alignments','250']
    result=subprocess.run(command,capture_output=True,text=True,timeout=120);(d/(method['name']+'.raw.tsv')).write_text(result.stdout);(d/(method['name']+'.stderr')).write_text(result.stderr);assert result.returncode==0
    rows=[dict(zip(m.FIELDS,x.split('\t'))) for x in result.stdout.splitlines()]
    if translated_query:rows=[t.normalize(x,mapping,records) for x in rows]
    summary={}
    for q in records:
        hits=[x for x in rows if x['qseqid']==q];intended=[x for x in hits if x['qframe']=='1' and x['sframe']=='1'] if method['executable']!='blastn' else hits
        summary[q]=dict(any_hit=bool(hits),intended_frame_hit=bool(intended),best_evalue=min((float(x['evalue']) for x in intended),default=None),best_bitscore=max((float(x['bitscore']) for x in intended),default=None),max_coordinate_query_coverage_pct=max((100*(abs(int(x['qend'])-int(x['qstart']))+1)/len(records[q]) for x in intended),default=None))
    results[method['name']]=dict(command=command,returncode=result.returncode,rows=len(rows),queries=summary)
m.save(d/'characterization.json',dict(created_utc=m.now(),source_sha256=m.sha(__file__),reference_sha256=m.sha(f/'reference.fna'),query_sha256=m.sha(d/'queries.fna'),truth=truth,methods=results,limits=['One deterministic synthetic80AAreference with controlled synonymous recoding and AA substitutions','Not an estimate of recall or accuracy on real samples','Tiny-database Evalues cannot establish significance against comprehensive core_nt','No speed claims from this fixture','Both even and random mutation placement fixed before evaluating methods']))
print(json.dumps({method:{q:r['max_coordinate_query_coverage_pct'] for q,r in v['queries'].items()} for method,v in results.items()}))

"""Characterize no-hit denominators, strict filtering, disagreement and artificial taxa."""
import csv, importlib.util, json
from pathlib import Path
spec=importlib.util.spec_from_file_location('analysis',Path(__file__).with_name('analyze.py'));a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a);m=a.m
root=m.ROOT/'cohort-storage/translated-core-nt/analysis-fixture';root.mkdir(exist_ok=False)
queries={name+'/1':'ACGT'*25 for name in ['control','disagreement','weak','none','artificial']}
(root/'queries.fna').write_text(''.join('>'+q+'\n'+seq+'\n' for q,seq in queries.items()))
a.write(root/'selection.tsv',[dict(seq_id=q.rsplit('/',1)[0],sample='synthetic',stratum='fixture') for q in queries],['seq_id','sample','stratum'])
m.save(root/'input-manifest.json',dict(query_sha256=m.sha(root/'queries.fna')))
for method in ['megablast','blastn','tblastx']:
    d=root/method;d.mkdir();rows=[]
    for q in queries:
        if q=='none/1' or (q=='disagreement/1' and method=='megablast'):continue
        taxid='81077' if q=='artificial/1' else '2' if q=='disagreement/1' and method=='blastn' else '10239'
        row=dict(qseqid=q,sseqid='subject',staxids=taxid,bitscore='100',evalue='0.1' if q=='weak/1' else '1e-20',length='30',pident='100',qcovs='90',qcovhsp='90',qlen='100',slen='100',qstart='1',qend='90',sstart='1',send='90',qframe='1',sframe='1',qseq='A'*30,sseq='A'*30)
        rows.append(row)
    with (d/'alignments.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=m.FIELDS,delimiter='\t');w.writerows(rows)
    m.save(d/'complete.json',{});m.save(d/'resources.json',dict(status='complete',wall_seconds=1));(d/'time.txt').write_text('User time (seconds): 1\nSystem time (seconds): 0\nMaximum resident set size (kbytes): 100\n')
a.main(root)
rows=list(csv.DictReader((root/'analysis/per-mate.tsv').open(),delimiter='\t'));assert len(rows)==30
for x in rows:
    if x['query']=='none/1':assert x['status']=='no_alignment'
    if x['query']=='weak/1':assert x['status']==('no_alignment' if x['threshold']=='strict_E1e-5' else 'viral_only')
    if x['query']=='artificial/1':assert x['pipeline_natural_lca']=='NA' and x['pipeline_combined_lca']=='81077'
s=json.loads((root/'analysis/summary.json').read_text());assert all(x['disagreeing_mates']==1 for x in s['disagreements'].values())
m.save(root/'check-complete.json',dict(invariants=['All5queries retained across3methods and2thresholds','No-hit query never dropped','Evalue strict filter removes weak evidence','Artificial-only evidence has no natural LCA','Each comparison identifies the one changed query'],excluded_from_real_data=True));print('ANALYSIS_INVARIANTS_VERIFIED')

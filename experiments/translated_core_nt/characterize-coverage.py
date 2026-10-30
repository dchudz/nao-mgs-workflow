"""Separate translated HSP coverage filtering from displayed coverage values."""
import importlib.util,json,subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
f=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=f/'coverage';d.mkdir(exist_ok=False);records=[]
for threshold in [0,30,35,40,80,99]:
    cmd=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE,'tblastx','-db',str(f/'db'),'-query',str(f/'queries.fna'),'-outfmt',m.OUTFMT,'-qcov_hsp_perc',str(threshold),'-evalue','10','-max_hsps','5','-num_alignments','250']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=60);assert result.returncode==0
    (d/(str(threshold)+'.tsv')).write_text(result.stdout);rows=[]
    for line in result.stdout.splitlines():
        x=dict(zip(m.FIELDS,line.split('\t')));span=100*(abs(int(x['qend'])-int(x['qstart']))+1)/int(x['qlen'])
        rows.append(dict(query=x['qseqid'],frame=int(x['qframe']),aa_length=int(x['length']),coordinate_coverage_pct=span,reported_qcovhsp=float(x['qcovhsp']),query_start=int(x['qstart']),query_end=int(x['qend'])))
        assert span+0.51>=threshold,'Returned HSP violates the requested coverage threshold'
    records.append(dict(threshold=threshold,command=cmd,rows=rows))
# Full coding-frame match survives99percent. Known short HSP moves at its coordinate-derived boundary.
assert all(any(x['query']==q and x['frame']==1 and x['aa_length']==80 for x in records[-1]['rows']) for q in ['identical','synonymous'])
for threshold,expected in [(30,True),(35,False)]:
    observed=any(x['query']=='identical' and x['aa_length']==26 for r in records if r['threshold']==threshold for x in r['rows']);assert observed==expected
m.save(d/'characterization.json',dict(created_utc=m.now(),records=records,result='Requested filter respects coordinate-derived coverage on these fixtures; reported qcovhsp is inconsistent when multiple translated HSPs are present. Preserve raw values and derive coverage from nucleotide coordinates for interpretation.',limits='Finite characterization, not proof of all BLAST edge cases'))
print(json.dumps({str(r['threshold']):r['rows'] for r in records}))

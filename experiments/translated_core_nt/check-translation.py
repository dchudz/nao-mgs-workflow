"""Executable characterization: synonymous coding preserves protein identity."""
import importlib.util, json, random, subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
d=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d.mkdir(exist_ok=False)
codons=[('GCT','GCC'),('CGT','CGC'),('AAT','AAC'),('GAT','GAC'),('TGT','TGC'),('CAA','CAG'),('GAA','GAG'),('GGT','GGC'),('CAT','CAC'),('ATT','ATC'),('CTT','CTC'),('AAA','AAG'),('TTT','TTC'),('CCT','CCC'),('TCT','TCC'),('ACT','ACC'),('TAT','TAC'),('GTT','GTC')]
rng=random.Random(20261022);pairs=[rng.choice(codons) for _ in range(80)];reference=''.join(p[0] for p in pairs);recoded=''.join(p[1] for p in pairs)
assert sum(a==b for a,b in zip(reference,recoded))==160
(d/'reference.fna').write_text('>reference\n'+reference+'\n');(d/'queries.fna').write_text('>identical\n'+reference+'\n>synonymous\n'+recoded+'\n')
prefix=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE]
subprocess.run(prefix+['makeblastdb','-in',str(d/'reference.fna'),'-dbtype','nucl','-out',str(d/'db')],check=True,stdout=subprocess.DEVNULL)
results={}
for name in ['megablast','blastn','tblastx']:
    config=json.loads(Path(__file__).with_name(name+'.json').read_text())
    cmd=prefix+[config['executable']]+config['arguments']+['-db',str(d/'db'),'-query',str(d/'queries.fna'),'-out',str(d/(name+'.tsv')),'-evalue','10','-qcov_hsp_perc','30','-max_hsps','5','-num_alignments','250','-num_threads','1','-outfmt',m.OUTFMT]
    subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
    rows=[dict(zip(m.FIELDS,line.split('\t'))) for line in (d/(name+'.tsv')).read_text().splitlines()]
    assert any(x['qseqid']=='identical' for x in rows)
    if name=='tblastx':assert any(x['qseqid']=='synonymous' and float(x['pident'])==100 and int(x['length'])==80 and x['qseq']==x['sseq'] for x in rows)
    results[name]=dict(aligned_queries=sorted({x['qseqid'] for x in rows}),rows=len(rows),command=cmd)
m.save(d/'complete.json',dict(completed_utc=m.now(),nucleotide_identity_pct=100*160/240,translated_identical_residues=80,results=results,excluded_from_real_data=True))
print(json.dumps({name:{k:v for k,v in value.items() if k!='command'} for name,value in results.items()}))

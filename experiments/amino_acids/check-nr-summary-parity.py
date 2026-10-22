"""Preserve the original near-best taxon evidence for unsorted, redundant alignments."""
import importlib.util,json,random
from pathlib import Path
r=Path(__file__).resolve().parent/'experiment-code'
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod.summarize
old=load('original',r/'protein-nr-summary.original.py');new=load('bounded',r/'protein-nr-summary.py')
parents={1:1,10239:1,10:10239,11:10,12:10,2:1,21:2};rng=random.Random(20261021)
rows=[];calls={f'p{i}':rng.choice([10,11,12,21,0]) for i in range(50)}
for i in range(10000):
 q=rng.choice(list(calls))+rng.choice(['/1','/2']);taxids=rng.choice(['11','12','11;12','21','11;21','0','999','11;999',''])
 rows.append([q,'accession',str(rng.choice([89.9,90,90.1,100,120,200])),'1e-8',str(rng.choice([24,25,50])),'90',taxids])
rows.extend([['missing/1','x','100','1e-8','24','90','11']]);calls['missing']=11;calls['absent']=11
for order in [rows,rows[::-1],rng.sample(rows,len(rows))]:assert old(iter(order),parents,calls)==new(iter(order),parents,calls)
print(json.dumps({'alignment_rows':len(rows),'query_pairs':len(calls),'orders_checked':3,'exact_result_parity':True,'coverage':['minimum alignment length','near-best boundary','repeated accession taxa','higher scores arriving later','mixed viral/nonviral','missing taxonomy','absent alignments','one/two mates']}))

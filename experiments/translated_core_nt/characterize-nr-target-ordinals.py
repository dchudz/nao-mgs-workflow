"""Check output-column parity and target ordinals on an isolated synthetic native BLAST database."""
import argparse, collections, hashlib, json, random, subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--commit',required=True);a=p.parse_args()
r=Path('/home/ec2-user/experiment');d=r/'protein/nr-target-ordinal-fixture';d.mkdir(exist_ok=False)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
exe=r/'protein/tools/diamond';rng=random.Random(17);alphabet='ACDEFGHIKLMNPQRSTVWY'
peptides={name:''.join(rng.choice(alphabet) for _ in range(80)) for name in ('CONTROL_A','CONTROL_B')}
codons=dict(zip(alphabet,[('GCT','GCC'),('TGT','TGC'),('GAT','GAC'),('GAA','GAG'),('TTT','TTC'),('GGT','GGC'),('CAT','CAC'),('ATT','ATC'),('AAA','AAG'),('CTG','CTC'),('ATG','ATG'),('AAT','AAC'),('CCT','CCC'),('CAA','CAG'),('CGT','CGC'),('TCT','TCC'),('ACT','ACC'),('GTT','GTC'),('TGG','TGG'),('TAT','TAC')]))
encode=lambda peptide,option:''.join(codons[c][option] for c in peptide)
dna=encode(peptides['CONTROL_A'],0)
queries={'A':dna,'A_synonymous':encode(peptides['CONTROL_A'],1),'A_reverse':dna.translate(str.maketrans('ACGT','TGCA'))[::-1],'B':encode(peptides['CONTROL_B'],0)}
for name,seqs in [('targets.faa',peptides),('queries.fna',queries)]:
    (d/name).write_text(''.join('>'+key+'\n'+value+'\n' for key,value in seqs.items()))
make=['docker','run','--rm','--network','none','--cpus','1','--memory','2g','-v',str(d)+':'+str(d),image,'makeblastdb','-in',str(d/'targets.faa'),'-dbtype','prot','-parse_seqids','-taxid','1','-out',str(d/'db')]
with (d/'makeblastdb.log').open('w') as log:subprocess.run(make,check=True,stdout=log,stderr=log)
fields='qseqid sseqid bitscore evalue length pident staxids'.split()
commands=[];checks=[]
for mode in ('0','1'):
    outputs=[]
    for extra in (False,True):
        name='mode'+mode+('-ordinal' if extra else '-original')
        cmd=[str(exe),'blastx','--db',str(d/'db'),'--query',str(d/'queries.fna'),'--out',str(d/(name+'.tsv')),'--sensitive','--threads','2','--block-size','2','--index-chunks','4','--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--algo',mode,'--outfmt','6']+fields+(['snum'] if extra else [])
        commands.append(cmd)
        with (d/(name+'.log')).open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
        rows=[line.split('\t') for line in (d/(name+'.tsv')).read_text().splitlines()];assert rows and {x[0] for x in rows}==queries.keys()
        assert all(len(x)==(8 if extra else 7) for x in rows)
        outputs.append(rows)
    assert collections.Counter(tuple(x) for x in outputs[0])==collections.Counter(tuple(x[:7]) for x in outputs[1])
    identifiers=sorted({(x[1],int(x[7])) for x in outputs[1]})
    assert identifiers==[('CONTROL_A',0),('CONTROL_B',1)],identifiers
    checks.append(dict(mode=mode,queries=len(queries),rows=len(outputs[0]),first_seven_columns_multiset_equal=True,observed_subject_ordinals=identifiers))
result=dict(source_commit=a.commit,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),diamond_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),BLAST_container=image,make_command=make,commands=commands,checks=checks,limits=['Synthetic mechanism fixture only, excluded from biological validation and real-data counts.','Uses a separate two-sequence native BLAST protein database; the comprehensive NR snapshot is unchanged.','Tests that adding snum preserves the seven original output columns in both modes and reports zero-based target ordinals in this fixture.','Full-reference raw-column parity remains measured by the controlled repeat, not inferred from this fixture alone.'])
(d/'complete.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(checks),flush=True)

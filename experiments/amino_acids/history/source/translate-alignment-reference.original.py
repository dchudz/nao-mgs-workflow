"""Six-frame protein alignment reference from the exact Bowtie2 masked FASTA."""
import gzip,re,json,time,hashlib,random
from pathlib import Path
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Data import CodonTable
r=Path('/home/ec2-user/experiment')
base=np.full(256,4,dtype=np.uint8)
for i,b in enumerate(b'TCAG'):base[b]=i
lut=np.full(125,ord('*'),dtype=np.uint8)
for codon,aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items():
 a,b,c=[int(base[ord(ch)]) for ch in codon];lut[a*25+b*5+c]=ord(aa)
complement=bytes.maketrans(b'ACGT',b'TGCA')
def peptides(dna,minimum):
 dna=dna.upper().replace('U','T').encode()
 for strand,seq in enumerate((dna,dna.translate(complement)[::-1])):
  codes=base[np.frombuffer(seq,dtype=np.uint8)]
  for frame in range(3):
   n=(len(seq)-frame)//3
   if n<=0:continue
   trip=codes[frame:frame+3*n].reshape(-1,3)
   aa=lut[trip[:,0]*25+trip[:,1]*5+trip[:,2]].tobytes()
   pos=0
   for part in aa.split(b'*'):
    if len(part)>=minimum:yield strand*3+frame,frame+pos*3,part.decode()
    pos+=len(part)+1

def oracle(dna,minimum):
 dna=dna.upper().replace('U','T')
 for strand,seq in enumerate((dna,str(Seq(dna).reverse_complement()))):
  for frame in range(3):
   parts=[]
   for i in range(frame,len(seq)-2,3):
    codon=seq[i:i+3]
    parts.append(str(Seq(codon).translate()) if set(codon)<=set('ACGT') else '*')
   for match in re.finditer('[ACDEFGHIKLMNPQRSTVWY]+',''.join(parts)):
    if len(match.group())>=minimum:yield strand*3+frame,frame+match.start()*3,match.group()
g=random.Random(20261021)
for dna in ['ATGGCTGAATTT','GCNGCNGCNGCN','uauAUGgggtaa']+[''.join(g.choices('ACGTNRY',k=g.randrange(5,180))) for _ in range(500)]:
 assert list(peptides(dna,2))==list(oracle(dna,2)),dna
src=r/'reference/results/virus-genomes-masked.fasta.gz';dst=r/'protein/viral-alignment.faa'
start=time.time();counts={'genomes':0,'proteins':0,'residues':0};h=hashlib.sha256()
with gzip.open(src,'rt') as f,dst.open('xb') as out:
 for title,dna in SimpleFastaParser(f):
  name=title.split()[0];assert '~' not in name;counts['genomes']+=1
  for frame,offset,protein in peptides(dna,25):
   raw=f'>{name}~{frame}~{offset}\n{protein}\n'.encode();out.write(raw);h.update(raw)
   counts['proteins']+=1;counts['residues']+=len(protein)
counts.update(seconds=time.time()-start,sha256=h.hexdigest(),minimum_length=25,genetic_code=1,ambiguity='break',source=str(src),purpose='Protein alignment reference uses the same DUST-masked genomes as Bowtie2; the amino-acid k-mer index separately uses human-masked genomes.',oracle_cases=503)
(r/'protein/alignment-translation.json').write_text(json.dumps(counts,indent=2));print(json.dumps(counts),flush=True)

"""Independently check real DIAMOND HSP coordinates on synthetic translated DNA."""
import csv,hashlib,json,random,subprocess,sys
from pathlib import Path
from diamond_geometry import Hsp,geometry
R=Path('/home/ec2-user/experiment');D=R/'integration/real-diamond';D.mkdir(parents=True,exist_ok=True)
# One unambiguous standard-code codon per residue, chosen independently of the adapter.
codons=dict(zip('ACDEFGHIKLMNPQRSTVWY',['GCT','TGT','GAT','GAA','TTT','GGT','CAT','ATT','AAA','CTG','ATG','AAT','CCT','CAA','CGT','TCT','ACT','GTT','TGG','TAT']))
inverse={v:k for k,v in codons.items()};rng=random.Random(920)
protein=''.join(rng.choices(list(codons),k=60));dna=''.join(codons[a] for a in protein)
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
translate=lambda s:''.join(inverse[s[i:i+3]] for i in range(0,len(s),3))
references={};queries={}
for frame in range(6):
    offset=30+frame%3;oriented='N'*offset+dna+'N'*(300-offset-len(dna))
    references[f'g{frame}']=oriented if frame<3 else rc(oriented)
for frame in range(1,4):
    q='N'*(frame-1)+dna+'N'*7
    queries[f'q{frame}']=q;queries[f'q{-frame}']=rc(q)
queries['insertion']=dna[:75]+'GCT'+dna[75:]
queries['deletion']=dna[:75]+dna[78:]
(D/'reference.faa').write_text(''.join(f'>g{f}~{f}~{30+f%3}\n{protein}\n' for f in range(6)))
(D/'queries.fna').write_text(''.join(f'>{k}\n{v}\n' for k,v in queries.items()))
exe=R/'protein/tools/diamond';commands=[
    [str(exe),'makedb','--in',str(D/'reference.faa'),'--db',str(D/'reference')],
    [str(exe),'blastx','--db',str(D/'reference'),'--query',str(D/'queries.fna'),'--out',str(D/'hits.tsv'),
     '--very-sensitive','--algo','0','--threads','4','--block-size','2','--index-chunks','4','--min-orf','1',
     '--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10','--outfmt','6',
     'qseqid','sseqid','bitscore','evalue','length','pident','qlen','qstart','qend','sstart','send','qframe','slen','qseq_gapped','sseq_gapped']]
for i,cmd in enumerate(commands):
    with (D/f'command-{i}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
observed=set();gaps=0
with (D/'hits.tsv').open() as f:
    for r in csv.reader(f,delimiter='\t'):
        q,s,bits,evalue,length,identity,qlen,qs,qe,ss,se,qframe,slen,qa,sa=r
        h=Hsp(s,300,int(slen),int(qlen),int(qframe),int(qs),int(qe),int(ss),int(se),qa,sa);g=geometry(h)
        fragment=references[g.genome][g.ref_start:g.ref_end_exclusive]
        if int(s.split('~')[-2])>=3:fragment=rc(fragment)
        assert translate(fragment)==sa.replace('-',''),r
        query=queries[q];lo,hi=sorted([int(qs),int(qe)]);fragment=query[lo-1:hi]
        if int(qframe)<0:fragment=rc(fragment)
        assert translate(fragment)==qa.replace('-',''),r
        assert g.ref_start_unclipped<=g.ref_start and g.ref_end_unclipped>=g.ref_end_exclusive-1
        observed.add((q,int(s.split('~')[-2])));gaps+=int('-' in qa+sa)
expected={(q,f) for q in queries for f in range(6)}
assert observed==expected,(expected-observed)
assert gaps>0,'The fixture must exercise actual gapped translated alignments'
receipt={'verified':True,'query_reference_cases':len(observed),'gapped_hsps':gaps,'commands':commands,
         'diamond_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
         'scope':'Synthetic software integration only; no biological finding or reduced production reference.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

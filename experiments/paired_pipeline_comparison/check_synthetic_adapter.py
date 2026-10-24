"""Known-origin control using only artificial sequences and fixture taxonomy.

prepare freezes expectations without importing the adapter. run consumes that
fixture unchanged and exercises real DIAMOND, both adapter entry points, and
the production coordinate-then-similarity duplicate executables.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

CODONS = dict(A=('GCT','GCC'), C=('TGT','TGC'), D=('GAT','GAC'),
    E=('GAA','GAG'), F=('TTT','TTC'), G=('GGT','GGC'), H=('CAT','CAC'),
    I=('ATT','ATC'), K=('AAA','AAG'), L=('CTT','CTC'), M=('ATG','ATG'),
    N=('AAT','AAC'), P=('CCT','CCC'), Q=('CAA','CAG'), R=('CGT','CGC'),
    S=('TCT','TCC'), T=('ACT','ACC'), V=('GTT','GTC'), W=('TGG','TGG'),
    Y=('TAT','TAC'))
# Literal genomic intervals, independently laid out on 900-base references.
# Per frame: origin A mate1/mate2, origin B mate1/mate2; zero-based half-open.
INTERVALS = [
    [(60,180),(210,330),(360,480),(450,570)],
    [(61,181),(211,331),(361,481),(451,571)],
    [(62,182),(212,332),(362,482),(452,572)],
    [(720,840),(570,690),(420,540),(330,450)],
    [(719,839),(569,689),(419,539),(329,449)],
    [(718,838),(568,688),(418,538),(328,448)],
]

def rc(s):
    return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2)+'\n')

def write_tsv(p, fields, rows):
    with (gzip.open(p,'wt') if str(p).endswith('.gz') else open(p,'w')) as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        w.writeheader();w.writerows(rows)

def read_tsv(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
        return list(csv.DictReader(f,delimiter='\t'))

def prepare(d):
    d.mkdir(parents=True,exist_ok=False)
    ref=d/'reference'/'results';ref.mkdir(parents=True)
    rng=random.Random(2026102401)
    expected={};fastq=[];genomes=[];peptides=[];metadata=[];taxa=[];lengths={}
    nodes=['1\t|\t1\t|\tno rank\t|','81077\t|\t1\t|\tno rank\t|']
    names=['1\t|\troot\t|','81077\t|\tartificial sequences\t|']
    for frame, bounds in enumerate(INTERVALS):
        gid=f'synthetic_{frame}';taxid=str(900000000+frame);phase=frame%3
        peptide=''.join(rng.choice(list(CODONS)) for _ in range(180))
        gene=''.join(CODONS[x][0] for x in peptide)
        offset=30+phase
        oriented='N'*offset+gene+'N'*(900-offset-len(gene))
        genome=oriented if frame<3 else rc(oriented)
        genomes.append(f'>{gid}\n{genome}\n')
        peptides.append(f'>{gid}~{frame}~{offset}\n{peptide}\n')
        lengths[gid]=900;metadata.append(dict(genome_id=gid,taxid=taxid,species_taxid=taxid))
        taxa.append(dict(taxid=taxid));nodes.append(f'{taxid}\t|\t1\t|\tspecies\t|')
        names.append(f'{taxid}\t|\tSynthetic control {frame}\t|')
        for origin, starts, intervals in [('a',(30,180),bounds[:2]),('b',(330,420),bounds[2:])]:
            for variant in (('high','low','synonymous') if origin=='a' else ('high',)):
                sid=f'{gid}_{origin}_{variant}';mate_seqs=[];mate_expect=[]
                for mate,(start,(lo,hi)) in enumerate(zip(starts,intervals)):
                    dna=gene[start:start+120]
                    # The literal expected interval must extract the independently
                    # generated coding sequence, before any adapter is imported.
                    assert (genome[lo:hi] if frame<3 else rc(genome[lo:hi]))==dna
                    if variant=='synonymous':
                        aa=peptide[start//3:(start+120)//3]
                        dna=''.join(CODONS[x][i%7==0] for i,x in enumerate(aa))
                        reverse_code={codon:letter for letter,cs in CODONS.items() for codon in cs}
                        assert ''.join(reverse_code[dna[i:i+3]] for i in range(0,120,3))==aa
                        assert dna!=gene[start:start+120]
                    seq='N'*phase+dna+'NNN'
                    if mate:seq=rc(seq)
                    qual=('I' if variant=='high' else ('5' if variant=='low' else '?'))*len(seq)
                    fastq.append(f'@{sid}/{mate+1}\n{seq}\n+\n{qual}\n')
                    left,right=(phase,3) if frame<3 else (3,phase)
                    mate_seqs.append(dict(sequence=seq,quality=qual))
                    mate_expect.append(dict(ref_start=lo,ref_end_exclusive=hi,
                        ref_start_unclipped=lo-left,ref_end_unclipped=hi-1+right,
                        query_rc=(frame>=3) != bool(mate)))
                expected[sid]=dict(genome=gid,taxid=taxid,origin=origin,
                    exemplar=f'{gid}_{origin}_high',mates=mate_expect,reads=mate_seqs)
    for mate in (1,2):fastq.append(f'@negative/{mate}\n'+ 'N'*123+'\n+\n'+'I'*123+'\n')
    (d/'sample.fastq').write_text(''.join(fastq))
    (d/'genomes.fna').write_text(''.join(genomes));(d/'peptides.faa').write_text(''.join(peptides))
    save(d/'genome-lengths.json',lengths)
    write_tsv(ref/'virus-genome-metadata-gid.tsv.gz',['genome_id','taxid','species_taxid'],metadata)
    write_tsv(ref/'total-virus-db-annotated.tsv.gz',['taxid'],taxa)
    (ref/'taxonomy-nodes.dmp').write_text('\n'.join(nodes)+'\n')
    (ref/'taxonomy-names.dmp').write_text('\n'.join(names)+'\n')
    save(d/'expected.json',dict(input_pairs=25,positive_pairs=24,exact_unique_pairs=19,
        alignment_exemplars=12,joint_exemplars=12,negative_ids=['negative'],rows=expected,
        coordinate_convention='Aligned intervals are zero-based half-open; unclipped end is zero-based inclusive.',
        scope='Artificial coding strings, no organism-derived sequences. Minimal controlled peptide database, not a test of full reference construction, environmental sensitivity, FASTP or validation.'))
    save(d/'manifest.json',{str(p.relative_to(d)):sha(p) for p in sorted(d.rglob('*')) if p.is_file()})

def check_rows(rows, expected, marked=False):
    actual={r['seq_id']:r for r in rows}
    assert len(actual)==len(rows)==expected['positive_pairs']
    assert set(actual)==set(expected['rows'])
    for sid,e in expected['rows'].items():
        row=actual[sid]
        assert row['prim_align_genome_id_all']==e['genome'],(sid,'genome')
        assert row['aligner_taxid_lca']==e['taxid'],(sid,'taxon')
        assert row['protein_aligned_mates']=='2',(sid,'mates')
        for i,suffix in enumerate(('','_rev')):
            for field in ('ref_start','ref_start_unclipped','ref_end_unclipped','query_rc'):
                assert row['prim_align_'+field+suffix]==str(e['mates'][i][field]),(sid,field,suffix,row['prim_align_'+field+suffix],e['mates'][i][field])
            assert row['query_seq'+suffix]==e['reads'][i]['sequence']
            assert row['query_qual'+suffix]==e['reads'][i]['quality']
        if marked:assert row['prim_align_dup_exemplar']==e['exemplar'],(sid,'duplicate group')
    if marked:
        joint={sid for sid,r in actual.items() if r['prim_align_dup_exemplar']==sid and r['sim_dup_exemplar']==sid}
        assert joint=={e['exemplar'] for e in expected['rows'].values()},joint
        assert len(joint)==expected['joint_exemplars']

def geometry_cases():
    from diamond_geometry import Hsp,geometry
    # Fixed answers for subject residues 11..50 of six 900-base references.
    spans=[(60,180),(61,181),(62,182),(720,840),(719,839),(718,838)]
    count=0
    for frame,(lo,hi) in enumerate(spans):
        for qframe,qstart,qend,phase in [(1,1,120,0),(2,2,121,1),(3,3,122,2),
                                      (-1,126,7,0),(-2,125,6,1),(-3,124,5,2)]:
            h=Hsp(f'synthetic_{frame}~{frame}~{30+frame%3}',900,180,126,
                  qframe,qstart,qend,11,50,'A'*40,'A'*40)
            g=geometry(h)
            left,right=(phase,6-phase) if frame<3 else (6-phase,phase)
            assert (g.ref_start,g.ref_end_exclusive,g.ref_start_unclipped,g.ref_end_unclipped,g.query_rc)==(
                lo,hi,lo-left,hi-1+right,(frame>=3)!=(qframe<0))
            count+=1
    return count

def run(a):
    d=a.fixture.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();commands=[]
    def command(args):
        args=list(map(str,args));t=time.monotonic()
        result=subprocess.run(args,capture_output=True,text=True)
        commands.append(dict(command=args,seconds=time.monotonic()-t,returncode=result.returncode))
        (out/f'command-{len(commands):02}.log').write_text(result.stdout+result.stderr)
        save(out/'commands.json',commands)
        assert result.returncode==0,(args,result.stderr[-3000:])
    manifest=json.loads((d/'manifest.json').read_text())
    for name,checksum in manifest.items():assert sha(d/name)==checksum,name
    expected=json.loads((d/'expected.json').read_text())
    src=a.source.resolve();scripts=src/'experiments/amino_acids';sys.path.insert(0,str(scripts))
    from diamond_adapter import HIT_FIELDS
    from protein_queries import prepare as queries
    from union_pairs import pairs
    n_geometry=geometry_cases()
    query_result=queries(d/'sample.fastq',out/'queries.fna',out/'queries.tsv.gz',prefix='control_')
    assert query_result['exact_unique_pairs']==expected['exact_unique_pairs']
    command([a.diamond,'makedb','--in',d/'peptides.faa','--db',out/'database'])
    command([a.diamond,'blastx','--db',out/'database','--query',out/'queries.fna','--out',out/'hits.tsv',
        '--very-sensitive','--algo','0','--threads','2','--block-size','2','--index-chunks','4',
        '--min-orf','1','--query-gencode','1','--evalue','1e-5','--query-cover','30','--top','10',
        '--outfmt','6',*HIT_FIELDS])
    hits=out/'hits.tsv';hits.write_text(''.join(sorted(hits.read_text().splitlines(keepends=True),key=lambda x:x.split('\t')[0])))
    common=['--hits',hits,'--members',out/'queries.fna.members.tsv.gz','--query-metadata',out/'queries.tsv.gz',
        '--reads',d/'sample.fastq','--sample','synthetic_control','--reference',d/'reference','--genome-lengths',d/'genome-lengths.json']
    for name,extra in [('diamond_adapter',[]),('parallel_diamond_adapter',['--workers','2'])]:
        command([sys.executable,scripts/(name+'.py'),*common,*extra,'--output',out/(name+'.tsv.gz')])
        check_rows(read_tsv(out/(name+'.tsv.gz')),expected)
    assert read_tsv(out/'diamond_adapter.tsv.gz')==read_tsv(out/'parallel_diamond_adapter.tsv.gz')
    command([sys.executable,scripts/'parallel_diamond_adapter.py',*common,'--workers','2','--extract-only','--output',out/'candidates.fastq.gz'])
    assert [raw for _,raw in pairs(out/'candidates.fastq.gz')]==[raw for key,raw in pairs(d/'sample.fastq') if key[0]!=b'@negative']
    command([a.mark_duplicates,'-i',out/'parallel_diamond_adapter.tsv.gz','-o',out/'alignment.tsv.gz','-m',out/'groups.tsv.gz','-d','1','-n','2'])
    rows=sorted(read_tsv(out/'alignment.tsv.gz'),key=lambda r:r['seq_id'])
    write_tsv(out/'alignment-sorted.tsv.gz',list(rows[0]),rows)
    command([a.mark_similarity,'-i',out/'alignment-sorted.tsv.gz','-o',out/'joint.tsv.gz'])
    rows=read_tsv(out/'joint.tsv.gz');check_rows(rows,expected,marked=True)
    # Prove the checker itself rejects the three errors this control targets.
    rejected=[]
    for field,value in [('prim_align_ref_start',str(int(rows[0]['prim_align_ref_start'])+1)),
                        ('aligner_taxid_lca','1'),('prim_align_dup_exemplar','wrong_group')]:
        changed=copy.deepcopy(rows);changed[0][field]=value
        try:check_rows(changed,expected,marked=True)
        except AssertionError:rejected.append(field)
        else:raise AssertionError('Checker accepted injected error: '+field)
    record=dict(passed=True,geometry_cases=n_geometry,positive_pairs=24,negative_pairs_excluded=1,
        exact_unique_pairs=query_result['exact_unique_pairs'],alignment_exemplars=len({r['prim_align_dup_exemplar'] for r in rows}),
        joint_exemplars=12,rejected_error_injections=rejected,seconds=time.monotonic()-start,
        expectation_sha256=sha(d/'expected.json'),fixture_manifest_sha256=sha(d/'manifest.json'),
        source_commit=a.source_commit,source_hashes={str(p.relative_to(src)):sha(p) for p in [
            scripts/'diamond_adapter.py',scripts/'diamond_geometry.py',scripts/'translated_coordinates.py',
            scripts/'parallel_diamond_adapter.py',scripts/'protein_queries.py',
            src/'modules/local/lcaTsv/resources/usr/bin/lca_tsv.py']},
        binary_hashes={str(p):sha(p) for p in (a.diamond,a.mark_duplicates,a.mark_similarity)},
        scope=expected['scope'])
    save(out/'result.json',record);print(json.dumps(record,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
    q=sub.add_parser('prepare');q.add_argument('fixture',type=Path)
    q=sub.add_parser('run')
    for flag in ('fixture','output','source','diamond','mark-duplicates','mark-similarity'):q.add_argument('--'+flag,type=Path,required=True)
    q.add_argument('--source-commit',required=True)
    args=p.parse_args()
    if args.mode=='prepare':prepare(args.fixture)
    else:run(args)

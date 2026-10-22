"""Test six-frame query-only translation, coordinate mapping, and tblastn tasks."""
import csv,importlib.util,json,subprocess
from pathlib import Path
from Bio.Seq import Seq
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def translate_queries(source,dest):
    records=m.fasta(source);mapping={}
    with Path(dest).open('x') as out:
        for q,dna in records.items():
            for sign,seq in [(1,dna),(-1,str(Seq(dna).reverse_complement()))]:
                for offset in range(3):
                    n=(len(seq)-offset)//3
                    if n<=0:continue
                    aa=str(Seq(seq[offset:offset+3*n]).translate(table=1));name=q+'~frame'+str(sign*(offset+1))
                    assert name not in mapping
                    mapping[name]=dict(query=q,frame=sign*(offset+1),nucleotide_length=len(dna),protein_length=n)
                    out.write('>'+name+'\n'+aa+'\n')
    assert all(sum(x['query']==q for x in mapping.values())==6 for q in records if len(records[q])>=5)
    return records,mapping

def normalize(row,mapping,records):
    meta=mapping[row['qseqid']];q=meta['query'];frame=meta['frame'];length=meta['nucleotide_length'];start=int(row['qstart']);end=int(row['qend']);assert start<=end
    left=abs(frame)-1+3*(start-1);right=abs(frame)-1+3*end
    source=records[q] if frame>0 else str(Seq(records[q]).reverse_complement())
    assert str(Seq(source[left:right]).translate(table=1)).upper()==row['qseq'].upper(),'Query coordinate translation mismatch'
    result=dict(row,qseqid=q,qlen=str(length),qframe=str(frame),qstart=str(left+1 if frame>0 else length-left),qend=str(right if frame>0 else length-right+1))
    result['qcovhsp']=str(100*(right-left)/length)
    # Aggregate per-subject coverage must be recomputed after collecting all frames/HSPs for real data.
    result['qcovs']='NA'
    return result

if __name__=='__main__':
    fixture=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture';d=fixture/'query-only';d.mkdir(exist_ok=False)
    records,mapping=translate_queries(fixture/'queries.fna',d/'queries.faa');m.save(d/'query-map.json',mapping)
    prefix=['docker','run','--rm','--network','none','-v',f'{m.ROOT}:{m.ROOT}',m.IMAGE]
    results={}
    for task in ['tblastn','tblastn-fast']:
        command=prefix+['tblastn','-task',task,'-ungapped','-comp_based_stats','0','-seg','12 2.2 2.5','-soft_masking','false','-db_gencode','1','-db',str(fixture/'db'),'-query',str(d/'queries.faa'),'-outfmt',m.OUTFMT,'-evalue','10','-max_hsps','5','-num_alignments','250']
        result=subprocess.run(command,capture_output=True,text=True,timeout=60);(d/(task+'.raw.tsv')).write_text(result.stdout)
        rows=[normalize(dict(zip(m.FIELDS,line.split('\t'))),mapping,records) for line in result.stdout.splitlines()] if result.returncode==0 else []
        exact=[x for x in rows if x['qseqid']=='synonymous' and float(x['pident'])==100 and int(x['length'])==80]
        results[task]=dict(command=command,returncode=result.returncode,stderr=result.stderr,rows=len(rows),synonymous80aa_exact_hit=bool(exact),exact_hits=exact)
        with (d/(task+'.normalized.tsv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=m.FIELDS,delimiter='\t');w.writerows(rows)
    m.save(d/'characterization.json',dict(created_utc=m.now(),query_sha256=m.sha(fixture/'queries.fna'),translated_query_sha256=m.sha(d/'queries.faa'),input_queries=len(records),translated_queries=len(mapping),results=results,limits=['Synthetic mechanism check only','No cost measurement','No claimed score or Evalue parity with tblastx','All six full frames retained including stops; no ORF or length filter','Real-data query coverage must be normalized in nucleotide coordinates']))
    assert all(x['returncode']==0 and x['synonymous80aa_exact_hit'] for x in results.values())
    print(json.dumps({task:{k:v for k,v in row.items() if k not in ['command','exact_hits']} for task,row in results.items()}))

"""Union two pair screens by identity, sequence, quality and multiplicity."""
import argparse,collections,gzip,json
from pathlib import Path

def pairs(path):
    with open(path,'rb') as check: compressed=check.read(2)==b'\x1f\x8b'
    with (gzip.open(path,'rb') if compressed else open(path,'rb')) as f:
        while True:
            a=[f.readline() for _ in range(4)]
            if not a[0]:
                if any(a):raise ValueError('truncated record')
                break
            b=[f.readline() for _ in range(4)]
            if not all(a+b):raise ValueError('truncated pair')
            ids=[]
            for mate,lines in enumerate((a,b),1):
                if not lines[0].startswith(b'@') or not lines[2].startswith(b'+'):raise ValueError('invalid FASTQ')
                token=lines[0].split()[0]
                if token.endswith((b'/1',b'/2')):token=token[:-2]
                ids.append(token)
                if len(lines[1].rstrip())!=len(lines[3].rstrip()):raise ValueError('sequence/quality mismatch')
            if ids[0]!=ids[1]:raise ValueError('unpaired IDs')
            key=(ids[0],a[1].rstrip(),a[3].rstrip(),b[1].rstrip(),b[3].rstrip())
            yield key,b''.join(a+b)

def union(base,amino,output,rescue):
    baseline=collections.Counter();stats=collections.Counter()
    with open(output,'xb') as out,open(rescue,'xb') as extra:
        for key,raw in pairs(base):baseline[key]+=1;out.write(raw);stats['nucleotide_pairs']+=1
        for key,raw in pairs(amino):
            stats['amino_pairs']+=1
            if baseline[key]>0:baseline[key]-=1;stats['shared_pairs']+=1
            else:out.write(raw);extra.write(raw);stats['rescued_pairs']+=1
    stats['union_pairs']=stats['nucleotide_pairs']+stats['rescued_pairs']
    assert stats['amino_pairs']==stats['shared_pairs']+stats['rescued_pairs']
    return dict(stats)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('base');p.add_argument('amino');p.add_argument('output');p.add_argument('rescue');p.add_argument('stats')
    a=p.parse_args();s=union(a.base,a.amino,a.output,a.rescue);Path(a.stats).write_text(json.dumps(s,indent=2));print(json.dumps(s))

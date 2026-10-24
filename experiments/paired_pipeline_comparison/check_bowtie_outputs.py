"""Read-only pair-conservation audit of completed Bowtie2 task outputs."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import time
from pathlib import Path

def identities(path):
    counts=collections.Counter()
    with gzip.open(path,'rb') as f:
        while True:
            first=[f.readline() for _ in range(4)]
            if not first[0]:
                if any(first):raise ValueError('Truncated FASTQ: '+str(path))
                break
            second=[f.readline() for _ in range(4)]
            if not all(first+second):raise ValueError('Truncated pair: '+str(path))
            ids=[];fields=[]
            for lines in (first,second):
                if not lines[0].startswith(b'@') or not lines[2].startswith(b'+'):raise ValueError('Invalid FASTQ')
                token=lines[0].split()[0]
                if token.endswith((b'/1',b'/2')):token=token[:-2]
                ids.append(token)
                seq,qual=lines[1].rstrip(b'\r\n'),lines[3].rstrip(b'\r\n')
                if len(seq)!=len(qual):raise ValueError('Sequence/quality length mismatch')
                fields.extend((seq,qual))
            if ids[0]!=ids[1]:raise ValueError('Pair identity mismatch')
            counts[(ids[0],hashlib.sha256(b'\0'.join(fields)).digest())]+=1
    return counts

p=argparse.ArgumentParser();p.add_argument('comparison',type=Path);a=p.parse_args()
started=time.monotonic();records=[]
for stage in ('baseline-run','protein-run'):
    receipt=json.loads((a.comparison/(stage+'.json')).read_text());assert receipt['exit_code']==0
    with (a.comparison/stage/'trace.tsv').open() as f:tasks=list(csv.DictReader(f,delimiter='\t'))
    for task in tasks:
        work=Path(task['workdir']);inputs=list(work.glob('*_bowtie2_in.fastq.gz'))
        if not inputs:continue
        assert task['status']=='COMPLETED' and len(inputs)==1
        base=inputs[0].name.removesuffix('_in.fastq.gz')
        mapped=work/(base+'_mapped.fastq.gz');unmapped=work/(base+'_unmapped.fastq.gz')
        source=identities(inputs[0]);kept=identities(mapped);removed=identities(unmapped)
        if source!=kept+removed:
            raise ValueError(f'Pair conservation failed in {work}: missing={sum((source-kept-removed).values())}, extra={sum(((kept+removed)-source).values())}')
        sam=work/(base+'_mapped.sam.gz')
        with gzip.open(sam,'rb') as f:
            while f.read(8*1024**2):pass
        records.append({'stage':stage,'process':task['name'],'workdir':str(work),
            'input_pairs':sum(source.values()),'mapped_pairs':sum(kept.values()),'unmapped_pairs':sum(removed.values()),
            'id_sequence_quality_multiset_equal':True,'mapped_sam_gzip_crc_valid':True})
assert len(records)>=5,'Expected viral/human/other baseline tasks and human/other protein tasks'
record={'passed':True,'seconds':time.monotonic()-started,'tasks':records,
        'scope':'Pair IDs, mate sequences, qualities and multiplicity are conserved across the mapped/unmapped split; valid gzip EOF also checked. This does not validate alignment correctness or prove writer-error propagation is generally safe.'}
(a.comparison/'bowtie-output-integrity.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))

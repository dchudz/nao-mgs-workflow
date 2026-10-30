"""Compare saved legacy and bounded adapters on identical real DIAMOND HSPs."""
import csv,gzip,hashlib,json
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/adapter-bounded'
work=R/'integration/protein-workflow/protein/work'
def compare(before,after):
    h=hashlib.sha256();rows=0
    with gzip.open(before,'rb') as a,gzip.open(after,'rb') as b:
        while True:
            x,y=a.read(1024*1024),b.read(1024*1024)
            assert x==y,(before,after)
            if not x:break
            h.update(x);rows+=x.count(b'\n')
    return {'bytes_identical':True,'decompressed_sha256':h.hexdigest(),'lines':rows,'before':str(before),'after':str(after)}
extract=compare(next(work.glob('*/*/candidates.fastq.gz')),D/'candidates.fastq.gz')
final=compare(next(work.glob('*/*/unsorted.tsv.gz')),D/'final.unsorted.tsv.gz')
trace=list(csv.DictReader((work.parent/'trace.tsv').open(),delimiter='\t'))
receipt={'verified':True,'scope':'Existing real HSPs, same filtered original reads; no new database search.',
         'eligible_HSPs':1330176,'eligible_queries':378,'maximum_HSPs_one_query':24827,
         'extract':extract,'final':final,
         'original_trace':[r for r in trace if r['name'].split()[0] in ('PROTEIN_CANDIDATES','PROTEIN_ADAPTER')],
         'bounded_extract':json.loads((D/'extract-timing.json').read_text()),
         'bounded_final':json.loads((D/'final-timing.json').read_text()),
         'timing_limit':'Legacy ran in the pinned Python container; bounded timing ran with worker Python3.12. Timing is diagnostic, not a controlled pipeline speed comparison.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

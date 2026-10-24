"""Measure the frozen adapter on bounded real-output slices; do not change results."""
import csv,gzip,hashlib,json,sys,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');E=R/'experiments/amino_acids';sys.path.insert(0,str(E))
from diamond_adapter import read_hits
D=R/'integration/protein-parser-scale';D.mkdir(exist_ok=True)
source=R/'e2e/SRR37076045/protein-run/work/62/f6a71d0c492f2c3e0e47e2f547717e/protein-hits.tsv'
lengths=json.loads((R/'reference/genome-lengths.json').read_text())
with gzip.open(R/'reference/results/total-virus-db-annotated.tsv.gz','rt') as f:viral={r['taxid'] for r in csv.DictReader(f,delimiter='\t')}
with gzip.open(R/'reference/results/virus-genome-metadata-gid.tsv.gz','rt') as f:genomes={r['genome_id']:(int(r['taxid'] if r['taxid'] in viral or r['species_taxid'] not in viral else r['species_taxid']),lengths[r['genome_id']]) for r in csv.DictReader(f,delimiter='\t')}
records=[]
with source.open('rb') as f:
    size=source.stat().st_size
    for part,fraction in enumerate([0,.3,.6,.9]):
        offset=int(size*fraction);f.seek(offset)
        if offset:f.readline()
        start=f.tell();target=D/f'part-{part}.tsv';lines=[]
        for _ in range(100000):
            line=f.readline();assert line.endswith(b'\n');lines.append(line)
        data=b''.join(lines);target.write_bytes(data)
        begin=time.perf_counter();cpu=time.process_time();eligible=sum(1 for _ in read_hits(target,genomes))
        seconds=time.perf_counter()-begin;cpu=time.process_time()-cpu
        record={'part':part,'byte_start':start,'physical_rows':len(lines),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'eligible_HSPs':eligible,'wall_seconds':seconds,'cpu_seconds':cpu,'full_current_bytes_scaled_single_pass_hours':size/len(data)*seconds/3600}
        records.append(record);print(json.dumps(record),flush=True)
receipt={'source':str(source),'source_bytes_at_sampling':size,'frozen_adapter_commit':(R/'source-commit.txt').read_text().strip(),'interpreter':sys.version,'parts':records,'scope':'Diagnostic throughput only: four bounded slices of real output, no search, no change to filtering or geometry. Scaled time is conditional on representativeness; final adapter also computes taxonomy and per-query summaries.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')

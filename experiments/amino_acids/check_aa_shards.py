"""Compare ordered independent AA-screen shards with the frozen whole-prefix result."""
import concurrent.futures,hashlib,json,resource,shutil,subprocess,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');old=R/'integration/aa-input-throughput';D=R/'integration/aa-shards'
D.mkdir(exist_ok=True);binary=R/'rust-tools/target/release/screen_amino_acids';index=R/'reference/virus.aa12.bin'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
count=5_000_000;size=625_000;parts=8
started=time.monotonic()
for mate in (1,2):
    with (old/f'SRR37076045_{mate}.fastq').open('rb') as f:
        for part in range(parts):
            target=D/f'{part:02d}.R{mate}.fastq'
            with target.open('wb') as out:
                for i in range(size):
                    h=f.readline();seq=f.readline();plus=f.readline();qual=f.readline()
                    assert h.rstrip()==f'@SRR37076045.{part*size+i+1}/{mate}'.encode()
                    assert len(seq.rstrip())==len(qual.rstrip())>0 and plus.startswith(b'+')
                    out.writelines([h,seq,plus,qual])
        assert not f.read(1)
prepare_seconds=time.monotonic()-started

def run(part):
    output=D/f'{part:02d}.matched.fastq'
    assert not output.exists(),'Do not overwrite a completed trial'
    cmd=[str(binary),'screen','--index',str(index),'--reads',str(D/f'{part:02d}.R1.fastq'),
         '--reads2',str(D/f'{part:02d}.R2.fastq'),'--threads','4','--matched',str(output)]
    start=time.monotonic()
    with (D/f'{part:02d}.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    return {'part':part,'command':cmd,'seconds':time.monotonic()-start,'output_sha256':sha(output),
            'counts':(D/f'{part:02d}.log').read_text().strip()}
start=time.monotonic();cpu=resource.getrusage(resource.RUSAGE_CHILDREN)
with concurrent.futures.ThreadPoolExecutor(parts) as pool:results=list(pool.map(run,range(parts)))
seconds=time.monotonic()-start;endcpu=resource.getrusage(resource.RUSAGE_CHILDREN)
output=D/'combined.matched.fastq'
with output.open('wb') as out:
    for part in range(parts):
        with (D/f'{part:02d}.matched.fastq').open('rb') as f:shutil.copyfileobj(f,out)
actual=sha(output);expected=sha(old/'plain.matched.fastq');assert actual==expected
receipt={'verified':True,'sample':'SRR37076045','pairs':count,'parts':parts,'pairs_per_part':size,
 'threads_per_part':4,'binary_sha256':sha(binary),'index_sha256':sha(index),'prepare_seconds':prepare_seconds,
 'screen_seconds':seconds,'screen_CPU_seconds':endcpu.ru_utime+endcpu.ru_stime-cpu.ru_utime-cpu.ru_stime,
 'output_sha256':actual,'baseline_output_sha256':expected,'trials':results,
 'scope':'Fixed 5-million-pair prefix, ordered concatenation must equal frozen whole-prefix AA12 output byte for byte. Includes index loading per shard; no scientific parameters changed.'}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='trials'}))

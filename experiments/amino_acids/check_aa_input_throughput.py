"""Compare input transport for the unchanged AA12 executable on a fixed prefix."""
import hashlib,json,resource,subprocess,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/aa-input-throughput';D.mkdir(parents=True,exist_ok=True)
sample='SRR37076045';tool=R/'sratoolkit.3.4.1-ubuntu64/bin/fastq-dump'
command=[str(tool),'--split-files','--skip-technical','--minReadLen','1','--defline-seq','@$ac.$si/$ri','--defline-qual','+',
         '-N','1','-X','5000000','-O',str(D),str(R/'data'/sample/sample)]
with (D/'decode.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
for mate in (1,2):
    p=D/f'{sample}_{mate}.fastq'
    with (D/f'{sample}_{mate}.fastq.gz').open('wb') as out:
        subprocess.run(['pigz','-1','-p','8','-c',str(p)],stdout=out,check=True)
results=[]
for name,extension in [('gzip-1','.gz'),('plain',''),('gzip-2','.gz')]:
    output=D/(name+'.matched.fastq');before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
    cmd=[str(R/'rust-tools/target/release/screen_amino_acids'),'screen','--index',str(R/'reference/virus.aa12.bin'),
         '--reads',str(D/f'{sample}_1.fastq{extension}'),'--reads2',str(D/f'{sample}_2.fastq{extension}'),
         '--threads','32','--matched',str(output)]
    with (D/(name+'.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'name':name,'command':cmd,'seconds':time.monotonic()-start,'CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'matched_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'counts':(D/(name+'.log')).read_text().strip()}
    assert 'total_pairs=5000000' in result['counts'];results.append(result)
    (D/'results.json').write_text(json.dumps({'sample':sample,'prefix_pairs':5000000,'input_command':command,'trials':results},indent=2)+'\n')
    print(json.dumps(result),flush=True)
assert len({r['matched_sha256'] for r in results})==1
(D/'verification.json').write_text(json.dumps({'verified':True,'output_bytes_identical':True,'sample':sample,'prefix_pairs':5000000,'trials':results},indent=2)+'\n')

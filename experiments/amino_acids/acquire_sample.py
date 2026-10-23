"""Acquire the frozen public FASTQ pair; verify both original ENA checksums."""
import concurrent.futures,datetime,hashlib,json,subprocess,sys,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');sample=sys.argv[1]
selection=json.loads((R/'experiments/amino_acids/end-to-end-full-depth/sample-selection.json').read_text())
row=next(x for x in selection['initial_samples'] if x['run_accession']==sample)
D=R/'data'/sample;D.mkdir(parents=True,exist_ok=True)
def get(i):
    target=D/f'R{i+1}.fastq.gz';partial=D/f'R{i+1}.partial'
    expected=int(row['fastq_bytes'].split(';')[i]);digest=row['fastq_md5'].split(';')[i]
    start=time.monotonic();url='https://'+row['fastq_ftp'].split(';')[i]
    if not target.exists():
        with (D/f'download-{i+1}.log').open('a') as log:
            subprocess.run(['curl','-fL','--retry','3','--connect-timeout','30','--max-time','14400','--continue-at','-',url,'-o',str(partial)],stdout=log,stderr=subprocess.STDOUT,check=True)
        source=partial
    else:source=target
    assert source.stat().st_size==expected,(source,'size')
    h=hashlib.md5()
    with source.open('rb') as f:
        for block in iter(lambda:f.read(16*1024**2),b''):h.update(block)
    assert h.hexdigest()==digest,(source,'md5')
    if source==partial:partial.rename(target)
    item={'mate':i+1,'source':url,'bytes':expected,'md5':digest,'seconds':time.monotonic()-start}
    print(json.dumps(item),flush=True);return item
with concurrent.futures.ThreadPoolExecutor(2) as pool:files=list(pool.map(get,range(2)))
receipt={'sample':sample,'files':files,'metadata':row,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'observed_pair_count':'Pending production COUNT_READS; metadata is not treated as a measured count.'}
(D/'verified.json').write_text(json.dumps(receipt,indent=2)+'\n')
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(D/'verified.json'),f's3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/acquisition/{sample}.json','--only-show-errors'],check=True)

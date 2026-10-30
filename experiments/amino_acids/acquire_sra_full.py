"""Decode every spot from the public normalized full-quality SRA object.

Bounded contiguous chunks are execution units, not sampling. Each original spot
ID must occur exactly once in both final FASTQ mates. SRA Lite is not used.
"""
import concurrent.futures,datetime,gzip,hashlib,json,shutil,subprocess,sys,time,xml.etree.ElementTree as ET
from pathlib import Path
R=Path('/home/ec2-user/experiment');sample=sys.argv[1]
rows=json.loads((R/'experiments/amino_acids/end-to-end-full-depth/sample-selection.json').read_text())['initial_samples']
row=next(r for r in rows if r['run_accession']==sample);D=R/'data'/sample;D.mkdir(parents=True,exist_ok=True)
T=R/'sratoolkit.3.4.1-ubuntu64/bin';AWS=['aws','--profile','pipeline-codex','--region','us-east-1']
source=f's3://sra-pub-run-odp/sra/{sample}/{sample}';archive=D/sample
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
def run(cmd,log):
    start=time.monotonic()
    with Path(log).open('w') as f:subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=14400)
    return time.monotonic()-start
meta=json.loads(subprocess.check_output(['aws','s3api','head-object','--no-sign-request','--bucket','sra-pub-run-odp','--key',f'sra/{sample}/{sample}'],text=True))
if not archive.exists():
    partial=D/(sample+'.partial')
    run(['aws','s3','cp',source,partial,'--no-sign-request','--only-show-errors'],D/'sra-download.log')
    assert partial.stat().st_size==meta['ContentLength'];partial.rename(archive)
assert archive.stat().st_size==meta['ContentLength']
if not (D/'sra-validated.json').exists():
    elapsed=run([T/'vdb-validate',archive],D/'vdb-validate.log')
    info={'public_object':source,'source_metadata':meta,'sha256':digest(archive),'vdb_validate_seconds':elapsed,'full_quality':True}
    (D/'sra-validated.json').write_text(json.dumps(info,indent=2)+'\n')
run([T/'sra-stat','--quick','--xml',archive],D/'sra-stat.xml')
xml=ET.parse(D/'sra-stat.xml');spots=int(xml.getroot().attrib['spot_count'])
assert spots==int(row['read_count']),(spots,row['read_count'])
chunk_size=5_000_000;chunks=D/'chunks';chunks.mkdir(exist_ok=True)
def decode(bounds):
    lo,hi=bounds;target=chunks/f'{lo:012d}';target.mkdir(exist_ok=True)
    marker=target/'complete.json'
    if marker.exists():return json.loads(marker.read_text())
    # Only partial generated chunk files are removed for a failed attempt.
    for p in target.glob('*.fastq*'):p.unlink()
    cmd=[T/'fastq-dump','--split-files','--skip-technical','--minReadLen','1',
         '--defline-seq','@$ac.$si/$ri','--defline-qual','+','-N',str(lo),'-X',str(hi),'-O',target,archive]
    elapsed=run(cmd,target/'decode.log');counts=[]
    for mate in (1,2):
        path=target/f'{sample}_{mate}.fastq';count=0
        # Verify every ID, quality length, and pair position before compression.
        with path.open('rb') as f:
            while True:
                h=f.readline()
                if not h:break
                seq,plus,qual=f.readline(),f.readline(),f.readline()
                assert h.rstrip()==f'@{sample}.{lo+count}/{mate}'.encode(),h
                assert plus.startswith(b'+') and len(seq.rstrip())==len(qual.rstrip())>0
                count+=1
        assert count==hi-lo+1,(lo,hi,mate,count);counts.append(count)
        run(['pigz','-1','-p','2',path],target/f'compress-{mate}.log')
    result={'lo':lo,'hi':hi,'mate_counts':counts,'decode_seconds':elapsed}
    marker.write_text(json.dumps(result));print(json.dumps({'sample':sample,'chunk_complete':result}),flush=True);return result
ranges=[(lo,min(lo+chunk_size-1,spots)) for lo in range(1,spots+1,chunk_size)]
with concurrent.futures.ThreadPoolExecutor(8) as pool:results=list(pool.map(decode,ranges))
files=[]
for mate in (1,2):
    target=D/f'R{mate}.fastq.gz';partial=D/f'R{mate}.assembling'
    if not target.exists():
        with partial.open('wb') as out:
            for lo,hi in ranges:
                with (chunks/f'{lo:012d}'/f'{sample}_{mate}.fastq.gz').open('rb') as f:shutil.copyfileobj(f,out,16*1024**2)
        partial.rename(target)
    files.append({'mate':mate,'bytes':target.stat().st_size,'sha256':digest(target),'pairs':spots})
receipt={'sample':sample,'source':source,'input_type':'NCBI normalized SRA, original quality values; not SRA Lite',
         'source_metadata':meta,'spot_count':spots,'files':files,'chunks':results,
         'headers':'Canonical @accession.spot/mate, consistently shared by all arms.',
         'selection':'Every spot from1 through verified spot_count, inclusive; chunks are not downsampling.',
         'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(D/'verified.json').write_text(json.dumps(receipt,indent=2)+'\n')
subprocess.run(AWS+['s3','cp',str(D/'verified.json'),f's3://pipeline-lab-839993661825-us-east-1/runs/e2e-full-depth-20261023/acquisition/{sample}.json','--only-show-errors'],check=True)
# Completed generated chunks duplicate the verified assembled files; release scratch.
shutil.rmtree(chunks)

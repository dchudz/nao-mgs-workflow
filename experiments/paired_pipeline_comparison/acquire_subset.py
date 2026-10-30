"""Uniform nested sampling of paired records from a verified full-quality archive.

Chunking changes execution only. IDs are selected across the complete run before
decoding; no sequence or pipeline result participates in selection.
"""
import argparse, concurrent.futures, datetime, gzip, hashlib, json, random
import shutil, subprocess, time, xml.etree.ElementTree as ET
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--sample',required=True)
p.add_argument('--selection',type=Path,required=True);p.add_argument('--tools',type=Path,required=True)
p.add_argument('--selector',type=Path,required=True);a=p.parse_args()
selection=json.loads(a.selection.read_text());row=next(r for r in selection['selected'] if r['run_accession']==a.sample)
deep=a.sample in selection['depth_runs']; n=int(row['read_count']); small=5_000_000; large=25_000_000 if deep else small
D=a.root/'inputs'/a.sample;D.mkdir(parents=True,exist_ok=True);archive=D/a.sample
start=time.monotonic()
def run(cmd,log,timeout=7200):
    t=time.monotonic()
    with log.open('w') as f:subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=timeout)
    return time.monotonic()-t
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''):h.update(b)
    return h.hexdigest()
meta=json.loads(subprocess.check_output(['aws','s3api','head-object','--bucket','sra-pub-run-odp','--key',f'sra/{a.sample}/{a.sample}','--no-sign-request'],text=True))
assert archive.stat().st_size==meta['ContentLength']
validation_seconds=run([a.tools/'vdb-validate',archive],D/'vdb-validate.log')
run([a.tools/'sra-stat','--quick','--xml',archive],D/'sra-stat.xml')
assert int(ET.parse(D/'sra-stat.xml').getroot().attrib['spot_count'])==n
seed=hashlib.sha256((selection['seed']+':pairs:'+a.sample).encode()).hexdigest()
rng=random.Random(int(seed,16)); chosen=rng.sample(range(1,n+1),large)
shallow=set(rng.sample(chosen,small)) if deep else None
chosen.sort();chunks=D/'chunks';chunks.mkdir(exist_ok=True)
size=5_000_000; groups={}
for spot in chosen:
    lo=(spot-1)//size*size+1
    groups.setdefault(lo,[]).append((spot,1 if shallow is None or spot in shallow else 0))
del chosen,shallow
tasks=[]
for lo,items in groups.items():
    out=chunks/f'{lo:012d}';out.mkdir(exist_ok=True);ids=out/'ids.tsv'
    with ids.open('w') as f:
        for spot,flag in items:f.write(f'{spot}\t{flag}\n')
    tasks.append((lo,min(lo+size-1,n),len(items),out,sha(ids)))
del groups
def decode(task):
    lo,hi,k,out,ids_sha=task;marker=out/'complete.json'
    if marker.exists():
        result=json.loads(marker.read_text());assert result['ids_sha256']==ids_sha
        for f in result['files']:assert (out/f['name']).stat().st_size==f['bytes'] and sha(out/f['name'])==f['sha256']
        return result
    elapsed=run([a.tools/'fastq-dump','--split-files','--skip-technical','--minReadLen','1','--defline-seq','@$ac.$si/$ri','--defline-qual','+','-N',lo,'-X',hi,'-O',out,archive],out/'decode.log')
    cmd=[a.selector,out/'ids.tsv',out/f'{a.sample}_1.fastq',out/f'{a.sample}_2.fastq',out/'selected',a.sample,lo,hi,int(deep),k]
    counts=json.loads(subprocess.check_output(list(map(str,cmd)),text=True))
    files=[]
    for level in (['shallow','deep'] if deep else ['shallow']):
        for mate in [1,2]:
            q=out/f'selected.{level}_{mate}.fastq';run(['pigz','-1','-p','2','-f',q],out/f'{level}-{mate}-compress.log')
            z=q.with_suffix(q.suffix+'.gz');files.append({'name':z.name,'bytes':z.stat().st_size,'sha256':sha(z)})
    for mate in [1,2]:(out/f'{a.sample}_{mate}.fastq').unlink()
    result={'first':lo,'last':hi,'ids_sha256':ids_sha,'counts':counts,'decode_seconds':elapsed,'files':files}
    marker.write_text(json.dumps(result)+'\n');print(json.dumps({'chunk':lo,'counts':counts}),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(8) as pool:results=list(pool.map(decode,tasks))
assert sum(r['counts']['selected'] for r in results)==large
assert sum(r['counts']['shallow'] for r in results)==small
files=[]
for level,count in ([('shallow',small),('deep',large)] if deep else [('shallow',small)]):
    dest=D/level;dest.mkdir(exist_ok=True)
    for mate in [1,2]:
        target=dest/f'R{mate}.fastq.gz';partial=target.with_suffix('.partial')
        with partial.open('wb') as f:
            for lo,hi,k,out,ids_sha in tasks:
                with (out/f'selected.{level}_{mate}.fastq.gz').open('rb') as source:shutil.copyfileobj(source,f,16*1024**2)
        partial.rename(target);files.append({'level':level,'mate':mate,'path':str(target),'pairs':count,'bytes':target.stat().st_size,'sha256':sha(target)})
receipt={'sample':a.sample,'population_pairs':n,'seed_sha256':seed,'sampling':'Exact uniform random pairs without replacement; nested 5M within 25M when applicable','full_quality':True,'source_metadata':meta,'source_sha256':sha(archive),'vdb_validation_seconds':validation_seconds,'files':files,'chunks':results,'seconds_excluding_initial_download':time.monotonic()-start,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(D/'verified.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'verified':True,'sample':a.sample,'seconds':receipt['seconds_excluding_initial_download']}),flush=True)

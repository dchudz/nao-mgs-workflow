"""Archive the exact core_nt files after all measured searches, without modifying shared backups."""
import concurrent.futures,datetime,hashlib,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/translated-core-nt';source=r/'scratch/blast_db_190b4ae1';out=d/'reference-archive';out.mkdir(exist_ok=True)
assert (d/'protein232/tblastx-word4/complete.json').exists() or (d/'protein232/tblastx-word4/resources.json').exists(),'Do not compete with an unfinished measured search'
resource=json.loads((d/'protein232/tblastx-word4/resources.json').read_text());assert resource['status'] in ['complete','time_limit','failed']
a=['aws','--profile','pipeline-codex','--region','us-east-1'];bucket='pipeline-lab-839993661825-us-east-1';prefix='runs/translated-core-nt-20261022/frozen-core-nt/'
assert subprocess.check_output(a+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
subprocess.run(a+['s3api','head-bucket','--bucket',bucket,'--expected-bucket-owner','839993661825'],check=True,stdout=subprocess.DEVNULL)
objects=json.loads((d/'pilot24/reference-inventory.json').read_text());assert len(objects)==836 and sum(x['bytes'] for x in objects)==318551974004
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def archive(o):
 name=o['name'];assert '/' not in name;path=source/name;assert not path.is_symlink();stat=path.stat();assert stat.st_size==o['bytes'] and stat.st_mtime_ns==o['mtime_ns'];sha=digest(path)
 if 'sha256' in o:assert sha==o['sha256']
 key=prefix+name;headcmd=a+['s3api','head-object','--bucket',bucket,'--key',key,'--expected-bucket-owner','839993661825','--checksum-mode','ENABLED'];existing=subprocess.run(headcmd,capture_output=True,text=True)
 if existing.returncode==0:target=json.loads(existing.stdout);uploaded=False
 else:
  assert '(404)' in existing.stderr or 'Not Found' in existing.stderr,existing.stderr
  subprocess.run(a+['s3','cp',str(path),'s3://'+bucket+'/'+key,'--storage-class','DEEP_ARCHIVE','--checksum-algorithm','CRC64NVME','--metadata','source-sha256='+sha,'--only-show-errors'],check=True);target=json.loads(subprocess.check_output(headcmd,text=True));uploaded=True
 assert path.stat().st_size==stat.st_size and path.stat().st_mtime_ns==stat.st_mtime_ns
 assert target['ContentLength']==o['bytes'] and target['StorageClass']=='DEEP_ARCHIVE' and target['ChecksumType']=='FULL_OBJECT' and target['ChecksumCRC64NVME'] and target['Metadata']['source-sha256']==sha
 return {'name':name,'bytes':o['bytes'],'source_sha256':sha,'destination_full_object_CRC64NVME':target['ChecksumCRC64NVME'],'uploaded_this_run':uploaded}
plan={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_script_sha256':digest(Path(__file__)),'original_inventory_sha256':digest(d/'pilot24/reference-inventory.json'),'source_files':len(objects),'source_bytes':sum(x['bytes'] for x in objects),'prefix':prefix,'monthly_data_storage_USD':sum(x['bytes'] for x in objects)/2**30*.00099,'minimum180day_data_storage_USD':sum(x['bytes'] for x in objects)/2**30*.00099*6,'verification':'Every source file SHA256 is recorded. AWS CLI uploads use end-to-end CRC64NVME checksums; stored full-object CRC and source-SHA metadata are read back. No archived payload is restored or independently downloaded for SHA verification.','limits':['Existing shared reference backup is not modified.','File sizes and modification times must match the frozen experiment inventory.','Archive restoration is required before reference reuse.','The failed or completed search status is preserved, never converted into a biological no-hit result.']};(out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan),flush=True)
records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(archive,o) for o in objects]
 try:
  for future in concurrent.futures.as_completed(futures):
   row=future.result();records.append(row)
   with (out/'verified-objects.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
   if len(records)%32==0:print(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archived_files':len(records),'total_files':len(objects),'archived_bytes':sum(x['bytes'] for x in records)}),flush=True)
 except BaseException:
  for future in futures:future.cancel()
  raise
assert len(records)==len(objects) and sum(x['bytes'] for x in records)==plan['source_bytes']
receipt=dict(plan,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),verified=True,files=sorted(records,key=lambda x:x['name']));(out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n')
for path in [out/'plan.json',out/'complete.json',d/'pilot24/reference-inventory.json',d/'pilot24/reference-info.txt']:
 target='s3://'+bucket+'/'+prefix+'receipts/'+path.name;subprocess.run(a+['s3','cp',str(path),target,'--only-show-errors'],check=True);download=subprocess.check_output(a+['s3','cp',target,'-','--only-show-errors']);assert hashlib.sha256(download).hexdigest()==digest(path)
print('FROZEN_CORE_NT_PRIVATE_ARCHIVE_CHECKSUM_VERIFIED',flush=True)

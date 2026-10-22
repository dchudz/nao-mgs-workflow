"""Retain the exact used NR snapshot privately after search completion, with whole-object CRC verification."""
import concurrent.futures,datetime,hashlib,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment/protein');root=r/'current-nr';manifest=json.loads((root/'source-manifest.json').read_text());out=r/'reference-archive';out.mkdir(exist_ok=True)
assert (r/'seed-nr-standalone/complete.json').exists(),'Do not compete with the measured NR search'
a=['aws','--profile','pipeline-codex','--region','us-east-1'];public=['aws','--no-sign-request','--region','us-east-1'];bucket='pipeline-lab-839993661825-us-east-1';prefix='runs/amino-acid-20261021/frozen-references/nr-20261017/'
assert subprocess.check_output(a+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
subprocess.run(a+['s3api','head-bucket','--bucket',bucket,'--expected-bucket-owner','839993661825'],check=True,stdout=subprocess.DEVNULL)
objects=manifest['objects'];assert len(objects)==1619 and sum(x['Size'] for x in objects)==918566558190
pricing=json.loads((r/'deep-archive-pricing.json').read_text());rate=next(float(x['rate']['pricePerUnit']['USD']) for x in pricing['rates'] if x['attributes']['usagetype']=='TimedStorage-GDA-ByteHrs');assert rate==.00099
plan={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_manifest_sha256':hashlib.sha256((root/'source-manifest.json').read_bytes()).hexdigest(),'source_bytes':manifest['bytes'],'objects':len(objects),'destination_bucket':bucket,'prefix':prefix,'storage_class':'DEEP_ARCHIVE','monthly_data_storage_USD':manifest['bytes']/2**30*rate,'minimum180day_data_storage_USD':manifest['bytes']/2**30*rate*6,'price_source':pricing['source'],'limits':['Private retention of already-used reference only, no new biological search.','Restore requires an archive restoration request; this is not an immediately readable reference copy.','Full-object CRC64NVME compares local uploads to public source objects whose ETags and sizes match the frozen manifest. It is not download-SHA verification.','Restore and request charges are separate from storage. No worker is stopped or deleted by this script.']};(out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan),flush=True)
def head(args,b,key):return json.loads(subprocess.check_output(args+['s3api','head-object','--bucket',b,'--key',key,'--checksum-mode','ENABLED'],text=True))
def upload(o):
 name=o['Key'][len(manifest['prefix']):];assert '/' not in name;path=root/name;assert path.is_file() and not path.is_symlink() and path.stat().st_size==o['Size'];stat=path.stat()
 source=head(public,manifest['bucket'],o['Key']);assert source['ETag']==o['ETag'] and source['ContentLength']==o['Size'] and source['ChecksumType']=='FULL_OBJECT';crc=source['ChecksumCRC64NVME'];key=prefix+name
 existing=subprocess.run(a+['s3api','head-object','--bucket',bucket,'--key',key,'--checksum-mode','ENABLED'],capture_output=True,text=True)
 if existing.returncode==0:target=json.loads(existing.stdout);uploaded=False
 else:
  assert '(404)' in existing.stderr or 'Not Found' in existing.stderr,existing.stderr
  subprocess.run(a+['s3','cp',str(path),'s3://'+bucket+'/'+key,'--storage-class','DEEP_ARCHIVE','--checksum-algorithm','CRC64NVME','--only-show-errors'],check=True);target=head(a,bucket,key);uploaded=True
 assert path.stat().st_size==stat.st_size and path.stat().st_mtime_ns==stat.st_mtime_ns
 assert target['ContentLength']==o['Size'] and target['StorageClass']=='DEEP_ARCHIVE' and target['ChecksumType']=='FULL_OBJECT' and target['ChecksumCRC64NVME']==crc
 return {'name':name,'bytes':o['Size'],'source_ETag':o['ETag'],'source_full_object_CRC64NVME':crc,'destination_full_object_CRC64NVME':target['ChecksumCRC64NVME'],'destination_ETag':target['ETag'],'uploaded_this_run':uploaded}
records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(upload,o) for o in objects]
 try:
  for future in concurrent.futures.as_completed(futures):
   row=future.result();records.append(row)
   with (out/'verified-objects.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
   if len(records)%32==0 or row['bytes']>10*2**30:print(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified_objects':len(records),'total_objects':len(objects),'verified_bytes':sum(x['bytes'] for x in records)}),flush=True)
 except BaseException:
  for future in futures:future.cancel()
  raise
assert len(records)==len(objects) and sum(x['bytes'] for x in records)==manifest['bytes']
receipt=dict(plan,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),verified=True,files=sorted(records,key=lambda x:x['name']))
(out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n')
for path in [out/'plan.json',out/'complete.json',root/'source-manifest.json',r/'deep-archive-pricing.json']:
 target='s3://'+bucket+'/'+prefix+'receipts/'+path.name;subprocess.run(a+['s3','cp',str(path),target,'--only-show-errors'],check=True);data=subprocess.check_output(a+['s3','cp',target,'-','--only-show-errors']);assert hashlib.sha256(data).digest()==hashlib.sha256(path.read_bytes()).digest()
print('FROZEN_NR_PRIVATE_ARCHIVE_CHECKSUM_VERIFIED',flush=True)

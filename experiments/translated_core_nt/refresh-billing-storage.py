"""Inventory current private experiment object classes without changing any object."""
import collections,datetime,json,os,subprocess
from pathlib import Path
p=Path(__file__).resolve().parent
if os.environ.get('PIPELINE_CLOSURE_USE_INHERITED_AWS_CREDENTIALS')=='1':
 env=os.environ.copy()
else:
 creds=json.loads(subprocess.check_output(['ssh','-F',str(p.parent/'aa-experiment-20261021/ssh_config'),'aa-experiment','aws configure export-credentials --profile pipeline-codex --format process'],text=True))
 env={k:v for k,v in os.environ.items() if not k.startswith('AWS_')}
 env.update(AWS_ACCESS_KEY_ID=creds['AccessKeyId'],AWS_SECRET_ACCESS_KEY=creds['SecretAccessKey'],AWS_SESSION_TOKEN=creds['SessionToken'])
env.update(AWS_CONFIG_FILE='/dev/null',AWS_SHARED_CREDENTIALS_FILE='/dev/null',AWS_EC2_METADATA_DISABLED='true',AWS_DEFAULT_REGION='us-east-1')
aws=['/opt/homebrew/bin/aws','--region','us-east-1']
def call(args):return json.loads(subprocess.check_output(aws+args,env=env,text=True))
assert call(['sts','get-caller-identity'])['Account']=='839993661825'
rates={'STANDARD':.023,'DEEP_ARCHIVE':.00099}
counts=collections.Counter();sizes=collections.Counter();prefixes=[]
for prefix in ['runs/amino-acid-20261021/','runs/translated-core-nt-20261022/']:
 rows=call(['s3api','list-objects-v2','--bucket','pipeline-lab-839993661825-us-east-1','--prefix',prefix,'--expected-bucket-owner','839993661825','--query','Contents[].[Size,StorageClass]'])
 pc=collections.Counter();pb=collections.Counter()
 for size,kind in rows:
  assert kind in rates,kind
  counts[kind]+=1;sizes[kind]+=size;pc[kind]+=1;pb[kind]+=size
 prefixes.append({'prefix':prefix,'objects':dict(pc),'bytes':dict(pb)})
record={'as_of_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'account':'839993661825','prefixes':prefixes,'objects_by_storage_class':dict(counts),'bytes_by_storage_class':dict(sizes),'USD_per_GiB_month':rates,'modeled_monthly_data_storage_USD':sum(sizes[k]/2**30*rates[k] for k in sizes),'deep_archive_180_day_minimum_data_storage_USD':sizes['DEEP_ARCHIVE']/2**30*rates['DEEP_ARCHIVE']*6,'pricing_source':'https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonS3/current/us-east-1/index.json','limits':['Current object versions only; old versions, incomplete multipart uploads, requests, transfer, retrieval and archive metadata overhead are excluded.','The exact two experiment prefixes exclude unrelated baseline runs.','Snapshot precedes the final small history/report backup; subsequent receipt/source files add a small amount.','Rates are the experiment pricing observations, not a new price quotation.']}
(p/'retained-storage-final.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)

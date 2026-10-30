"""Retire the exact owned primary only after reference/data/history verification."""
import datetime,hashlib,json,os,subprocess,time
from pathlib import Path
root=Path('/Users/david/Documents/ChatGPT/pipeline');d=root/'.tmp/translated-core-nt-20261022';work=Path('/Users/david/.codex/worktrees/amino-acid-screen/pipeline');ssh=root/'.tmp/aa-experiment-20261021/ssh_config'
account='839993661825';instance='i-061f366bd55def34c';volume='vol-098cdf9d650bd4a0d';bucket='pipeline-lab-839993661825-us-east-1'
closure_path=d/'primary-final-evidence/complete.json';closure=json.loads(closure_path.read_text());reference_path=d/'core-reference-archive-complete.json';reference=json.loads(reference_path.read_text());coverage=json.loads((d/'primary-symlink-coverage.json').read_text())
assert closure['account']==account and closure['all_scientific_archives_download_hash_verified']
assert {x['scope'] for x in closure['archives']}=={'cohort-storage','cohort','full','pilot','protein','reference-results','null-controls','metadata-and-inputs'}
assert reference['verified'] and len(reference['files'])==836 and reference['source_bytes']==318551974004
assert coverage['unpreserved_existing_scientific_targets']==[]
files={}
for name in closure['scope_receipts']:
 path=d/'primary-final-evidence'/name;record=json.loads(path.read_text());assert record['archive']['independent_download_hash_verified'];files.update({x['path']:x for x in record['files']})
required=['cohort-storage/protein-finalization/complete.json','cohort-storage/translated-core-nt/protein232/tblastx-word4/complete.json','cohort-storage/translated-core-nt/protein232/tblastx-word4-one-batch/resources.json','cohort-storage/translated-core-nt/protein232/analysis/summary.json','cohort-storage/translated-core-nt/protein232/pair-analysis/summary.json','cohort-storage/translated-core-nt/protein232/transition-analysis/summary.json']
for name in required:assert name in files,name
assert files['cohort-storage/translated-core-nt/protein232/tblastx-word4-one-batch/resources.json']['sha256']==hashlib.sha256((d/'protein232/tblastx-word4-one-batch/resources.json').read_bytes()).hexdigest()
head=subprocess.check_output(['git','-C',str(work),'rev-parse','HEAD'],text=True).strip()
assert not subprocess.check_output(['git','-C',str(work),'status','--porcelain'],text=True).strip()
backup=json.loads((d/'latest-source-backup-verification.json').read_text());assert backup['head']==head and backup['all_download_hashes_verified']
remote="python3 - <<'P'\nimport json,subprocess,hashlib,pathlib\nr=pathlib.Path('/home/ec2-user/experiment')\nassert not any(subprocess.run(['pgrep','-x',p],capture_output=True).stdout.strip() for p in ['blastn','tblastx','tblastn','diamond'])\np=r/'primary-closure-final/complete.json'\nprint(json.dumps({'complete_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'no_active_searches':True}))\nP"
remote_state=json.loads(subprocess.check_output(['ssh','-F',str(ssh),'aa-experiment',remote],text=True));assert remote_state['complete_sha256']==hashlib.sha256(closure_path.read_bytes()).hexdigest()
creds=json.loads(subprocess.check_output(['ssh','-F',str(ssh),'aa-experiment','aws configure export-credentials --profile pipeline-codex --format process'],text=True))
env={k:v for k,v in os.environ.items() if not k.startswith('AWS_')};env.update(AWS_ACCESS_KEY_ID=creds['AccessKeyId'],AWS_SECRET_ACCESS_KEY=creds['SecretAccessKey'],AWS_SESSION_TOKEN=creds['SessionToken'],AWS_CONFIG_FILE='/dev/null',AWS_SHARED_CREDENTIALS_FILE='/dev/null',AWS_EC2_METADATA_DISABLED='true',AWS_DEFAULT_REGION='us-east-1',PIPELINE_CLOSURE_USE_INHERITED_AWS_CREDENTIALS='1')
aws=['/opt/homebrew/bin/aws','--region','us-east-1']
def call(args):return json.loads(subprocess.check_output(aws+args,env=env,text=True))
assert call(['sts','get-caller-identity'])['Account']==account
block=call(['s3api','get-public-access-block','--bucket',bucket,'--expected-bucket-owner',account]);assert all(block['PublicAccessBlockConfiguration'].values())
base='s3://'+bucket+'/runs/translated-core-nt-20261022/'
for path,uri in [(closure_path,base+'primary-closure-final/complete.json'),(reference_path,'s3://'+bucket+'/'+reference['prefix']+'receipts/complete.json')]+[(d/'primary-final-evidence'/n,base+'primary-closure-final/'+n) for n in closure['scope_receipts']]:
 content=subprocess.check_output(aws+['s3','cp',uri,'-','--only-show-errors'],env=env);assert hashlib.sha256(content).hexdigest()==hashlib.sha256(path.read_bytes()).hexdigest()
for item in closure['archives']:
 assert item['uri'].startswith('s3://'+bucket+'/')
 h=call(['s3api','head-object','--bucket',bucket,'--key',item['uri'].removeprefix('s3://'+bucket+'/'),'--expected-bucket-owner',account]);assert h['ContentLength']==item['bytes']
for item in backup['files']:
 path=d/item['name'];assert path.stat().st_size==item['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256']
 h=call(['s3api','head-object','--bucket',bucket,'--key',(backup['prefix']+item['name']).removeprefix('s3://'+bucket+'/'),'--expected-bucket-owner',account]);assert h['ContentLength']==item['bytes']
state=call(['ec2','describe-instances','--instance-ids',instance])['Reservations'][0]['Instances'][0]
assert state['InstanceId']==instance and state['State']['Name']=='running' and state['InstanceType']=='x2iedn.4xlarge'
tags={x['Key']:x['Value'] for x in state['Tags']};assert tags['Name']=='pipeline-lab-aa-experiment' and tags['RunId']=='aa-20261021'
assert len(state['BlockDeviceMappings'])==1 and state['BlockDeviceMappings'][0]['Ebs']['VolumeId']==volume and state['BlockDeviceMappings'][0]['Ebs']['DeleteOnTermination']
receipt={'account':account,'instance':instance,'volume':volume,'requested_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_closure_sha256':hashlib.sha256(closure_path.read_bytes()).hexdigest(),'core_reference_objects':836,'source_backup_head':head,'symlink_coverage_sha256':hashlib.sha256((d/'primary-symlink-coverage.json').read_bytes()).hexdigest(),'public_access_block':block,'reason':'All bounded searches ended; exact reference, scientific archives and local experiment history are preserved and verified.'}
receipt['termination_response']=call(['ec2','terminate-instances','--instance-ids',instance]);out=d/'primary-worker-retirement.json';out.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'termination_requested_utc':receipt['requested_utc'],'instance':instance}),flush=True)
subprocess.run(aws+['ec2','wait','instance-terminated','--instance-ids',instance],env=env,check=True)
final=call(['ec2','describe-instances','--instance-ids',instance])['Reservations'][0]['Instances'][0];assert final['State']['Name']=='terminated';receipt['termination_confirmed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipt['final_state']=final['State']
for attempt in range(12):
 result=subprocess.run(aws+['ec2','describe-volumes','--volume-ids',volume],env=env,capture_output=True,text=True)
 if result.returncode and 'InvalidVolume.NotFound' in result.stderr:receipt['root_volume_deletion_confirmed']=True;break
 assert result.returncode==0,result.stderr
 time.sleep(5)
assert receipt.get('root_volume_deletion_confirmed')
other=call(['ec2','describe-instances','--filters','Name=tag:RunId,Values=aa-20261021']);receipt['owned_tagged_instance_states']=[{'id':i['InstanceId'],'state':i['State']['Name']} for r in other['Reservations'] for i in r['Instances']]
assert all(i['state']=='terminated' for i in receipt['owned_tagged_instance_states'])
out.write_text(json.dumps(receipt,indent=2)+'\n')
ret=d/'owned-worker-retirements.json';retired=json.loads(ret.read_text());assert instance not in retired;retired[instance]=receipt['termination_confirmed_utc'];ret.write_text(json.dumps(retired,indent=2)+'\n')
active=d/'active-runtime.json';s=json.loads(active.read_text());s.update(updated_utc=receipt['termination_confirmed_utc'],status='Experiments ended; both workers retired after verified private preservation',backup_refresh_needed=True);s['primary'].update(worker_state='terminated',root_volume_deleted=True,retirement_receipt=out.name,immediate_next_actions=[]);active.write_text(json.dumps(s,indent=2)+'\n')
# The inherited credentials remain in process memory for the final small receipt backup.
subprocess.run(['python3',str(d/'refresh-billing-storage.py')],env=env,check=True)
subprocess.run(['python3',str(root/'.tmp/aa-experiment-20261021/estimate-experiment-cost.py'),'--retirement-receipts',str(ret),'--output',str(d/'budget-ledger.json')],check=True)
ledger=json.loads((d/'budget-ledger.json').read_text());split=json.loads((d/'split-budget-ledger.json').read_text());split['as_of_utc']=ledger['as_of_utc']
for k,v in split['before_boundary'].items():split['after_boundary'][k]=ledger[k]-v
split['new_cap_modeled_remainder_excluding_other_charges']=500-split['after_boundary']['modeled_owned_compute_EBS_IPv4_total_USD'];(d/'split-budget-ledger.json').write_text(json.dumps(split,indent=2)+'\n')
subprocess.run(['python3',str(d/'write-experiment-results.py')],check=True)
results=work/'experiments/translated_core_nt/results/primary-retirement';results.mkdir(exist_ok=False)
for path in [out,ret,d/'split-budget-ledger.json',d/'budget-ledger.json',d/'retained-storage-final.json',d/'experiment-results.provenance.json']:(results/path.name).write_bytes(path.read_bytes())
subprocess.run(['git','-C',str(work),'add','experiments/translated_core_nt/results/primary-retirement'],check=True)
subprocess.run(['git','-C',str(work),'commit','-m','Record verified primary retirement and final split experiment spend'],check=True)
subprocess.run(['python3',str(d/'refresh-private-backup-local.py')],env=env,check=True)
print(json.dumps({'retired':instance,'root_volume_deleted':True,'all_tagged_workers_terminated':True,'second_window_modeled_USD':split['after_boundary']['modeled_owned_compute_EBS_IPv4_total_USD']}),flush=True)

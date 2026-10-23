"""Retire only the preserved NR experiment instance; retain auditable receipts."""
import datetime,hashlib,json,os,subprocess,time
from pathlib import Path
r=Path('/Users/david/Documents/ChatGPT/pipeline');d=r/'.tmp/translated-core-nt-20261022';ssh=r/'.tmp/aa-experiment-20261021/ssh_config'
account='839993661825';instance='i-0103f018916b6c9e5';volume='vol-07d4e85519abc9f76';bucket='pipeline-lab-839993661825-us-east-1'
closure_path=d/'nr-final-evidence-complete.json';closure=json.loads(closure_path.read_text());reference_path=d/'nr-reference-archive-complete.json';reference=json.loads(reference_path.read_text())
assert closure['account']==account and closure['archive_download_hash_verified'] and closure['frozen_NR_reference_objects']==1619
assert reference['verified'] and len(reference['files'])==1619 and reference['source_bytes']==918566558190
files={x['path']:x for x in closure['files']};raws={x['kind']:x for x in closure['previously_download_hash_verified_payloads']}
assert set(raws)=={'raw_NR','query_calls','native_controls','seeds'}
audit='protein/seed-nr-forced-double-index-trace/'
assert files[audit+'nr-with-ordinals.tsv']['sha256']=='0be98f28e7a03a3f517c892b77606e766a3dc93dff7dd1b90113c84641e3cfd9'
for name in ['output-audit-v2/complete.json','output-audit-v2/reference-label-inspection/complete.json','output-audit-v2/reference-label-inspection/sequence-reconstruction-verification.json','resources.json']:assert audit+name in files
symlinks=[]
for x in closure['excluded_symlinks']:
    assert x['target'].startswith('/home/ec2-user/experiment/')
    target=x['target'].removeprefix('/home/ec2-user/experiment/')
    covered=target in files or target=='protein/combined-validation/nr.tsv'
    assert covered,x
    symlinks.append(dict(x,covered_by='final archive' if target in files else raws['raw_NR']['uri']))
assert not subprocess.run(['ssh','-F',str(ssh),'aa-nr','pgrep -x diamond'],capture_output=True,text=True).stdout.strip()
creds=json.loads(subprocess.check_output(['ssh','-F',str(ssh),'aa-experiment','aws configure export-credentials --profile pipeline-codex --format process'],text=True))
env={k:v for k,v in os.environ.items() if not k.startswith('AWS_')};env.update(AWS_ACCESS_KEY_ID=creds['AccessKeyId'],AWS_SECRET_ACCESS_KEY=creds['SecretAccessKey'],AWS_SESSION_TOKEN=creds['SessionToken'],AWS_CONFIG_FILE='/dev/null',AWS_SHARED_CREDENTIALS_FILE='/dev/null',AWS_EC2_METADATA_DISABLED='true',AWS_DEFAULT_REGION='us-east-1')
aws=['/opt/homebrew/bin/aws','--region','us-east-1']
def call(args):return json.loads(subprocess.check_output(aws+args,env=env,text=True))
assert call(['sts','get-caller-identity'])['Account']==account
block=call(['s3api','get-public-access-block','--bucket',bucket,'--expected-bucket-owner',account]);assert all(block['PublicAccessBlockConfiguration'].values())
for path,uri in [(closure_path,'s3://'+bucket+'/runs/amino-acid-20261021/nr-closure-final/complete.json'),(reference_path,'s3://'+bucket+'/'+reference['prefix']+'receipts/complete.json')]:
    content=subprocess.check_output(aws+['s3','cp',uri,'-','--only-show-errors'],env=env);assert hashlib.sha256(content).hexdigest()==hashlib.sha256(path.read_bytes()).hexdigest()
for item in list(raws.values())+[dict(uri=closure['archive_uri'],bytes=closure['archive_bytes'])]:
    assert item['uri'].startswith('s3://'+bucket+'/')
    h=call(['s3api','head-object','--bucket',bucket,'--key',item['uri'].removeprefix('s3://'+bucket+'/'),'--expected-bucket-owner',account]);assert h['ContentLength']==item['bytes']
state=call(['ec2','describe-instances','--instance-ids',instance])['Reservations'][0]['Instances'][0]
assert state['InstanceId']==instance and state['State']['Name']=='running' and state['InstanceType']=='r7a.8xlarge'
tags={x['Key']:x['Value'] for x in state['Tags']};assert tags['Name']=='pipeline-lab-aa-nr' and tags['RunId']=='aa-20261021'
assert len(state['BlockDeviceMappings'])==1 and state['BlockDeviceMappings'][0]['Ebs']['VolumeId']==volume and state['BlockDeviceMappings'][0]['Ebs']['DeleteOnTermination']
receipt=dict(account=account,instance=instance,volume=volume,requested_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),final_archive_uri=closure['archive_uri'],final_archive_sha256=closure['archive_sha256'],frozen_reference_objects=1619,all_symlink_targets_preserved=symlinks,public_access_block=block,reason='All NR scientific outputs, raw data, exact reference, target inspection, failed attempts and source were preserved and verified. No NR search remains active.')
receipt['termination_response']=call(['ec2','terminate-instances','--instance-ids',instance])
(d/'nr-worker-retirement.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'instance':instance,'termination_requested_utc':receipt['requested_utc']}),flush=True)
subprocess.run(aws+['ec2','wait','instance-terminated','--instance-ids',instance],env=env,check=True)
final=call(['ec2','describe-instances','--instance-ids',instance])['Reservations'][0]['Instances'][0];assert final['State']['Name']=='terminated'
receipt['termination_confirmed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipt['final_state']=final['State']
for attempt in range(12):
    result=subprocess.run(aws+['ec2','describe-volumes','--volume-ids',volume],env=env,capture_output=True,text=True)
    if result.returncode and 'InvalidVolume.NotFound' in result.stderr:receipt['root_volume_deletion_confirmed']=True;break
    assert result.returncode==0,result.stderr
    time.sleep(5)
assert receipt.get('root_volume_deletion_confirmed'), 'Retired instance, but root-volume deletion needs inspection'
(d/'nr-worker-retirement.json').write_text(json.dumps(receipt,indent=2)+'\n')
path=d/'owned-worker-retirements.json';retired=json.loads(path.read_text()) if path.exists() else {};assert instance not in retired;retired[instance]=receipt['termination_confirmed_utc'];path.write_text(json.dumps(retired,indent=2)+'\n')
print(json.dumps({'instance':instance,'confirmed_utc':receipt['termination_confirmed_utc'],'root_volume_deletion_confirmed':True}),flush=True)

"""Recheck private selected-read objects whose local scratch copies were already released."""
import datetime,json,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');root=r/'cohort/heldout';prefix='runs/amino-acid-20261021/cohort-screens/';bucket='pipeline-lab-839993661825-us-east-1';expected={};records=list(root.glob('*/released-selected-scratch.json'))
for p in records:
 x=json.loads(p.read_text());base=prefix+p.parent.name+'/'
 assert x['private_s3_prefix']==f's3://{bucket}/'+base
 for name,size in x['verified_sizes'].items():
  assert '/' not in name and name.endswith('.fastq');expected[base+name]=size
aws=['aws','--profile','pipeline-codex','--region','us-east-1'];assert subprocess.check_output(aws+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
listing=json.loads(subprocess.check_output(aws+['s3api','list-objects-v2','--bucket',bucket,'--prefix',prefix,'--query','Contents[].{key:Key,bytes:Size}','--output','json'],text=True));actual={x['key']:x['bytes'] for x in listing}
missing=sorted(expected.keys()-actual.keys());mismatched=[k for k,v in expected.items() if k in actual and actual[k]!=v]
result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'release_records':len(records),'released_objects_checked':len(expected),'released_bytes_checked':sum(expected.values()),'missing':missing,'size_mismatches':mismatched,'verified':not missing and not mismatched,'scope':'Previously released local selected-read scratch objects remain privately backed up with their recorded byte sizes.'}
p=r/'screen-closure/released-backup-verification.json';p.write_text(json.dumps(result,indent=2));assert result['verified'],result
subprocess.run(aws+['s3','cp',str(p),f's3://{bucket}/runs/amino-acid-20261021/screen-closure/'+p.name,'--only-show-errors'],check=True);print(json.dumps(result))

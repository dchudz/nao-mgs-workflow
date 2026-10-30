"""Verify private S3 backup coverage and byte sizes before releasing a worker."""
import argparse,datetime,json,subprocess
from pathlib import Path

BUCKET='pipeline-lab-839993661825-us-east-1'
BASE='runs/amino-acid-20261021/'

def compare(expected,actual):
 missing=sorted(expected.keys()-actual.keys())
 mismatched=[{'key':k,'local_bytes':v,'s3_bytes':actual[k]} for k,v in expected.items() if k in actual and actual[k]!=v]
 return missing,mismatched

def main():
 p=argparse.ArgumentParser();p.add_argument('local_root');p.add_argument('prefix');p.add_argument('--exclude-dir',action='append',default=[]);p.add_argument('--output',required=True);a=p.parse_args()
 root=Path(a.local_root).resolve();prefix=a.prefix.rstrip('/')+'/'
 assert prefix.startswith(BASE) and prefix!=BASE and '..' not in prefix
 assert root.is_dir()
 aws=['aws','--profile','pipeline-codex','--region','us-east-1']
 account=subprocess.check_output(aws+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip();assert account=='839993661825'
 expected={};excluded=[]
 output=Path(a.output).resolve()
 for f in root.rglob('*'):
  rel=f.relative_to(root)
  if any(str(rel)==x or str(rel).startswith(x.rstrip('/')+'/') for x in a.exclude_dir):continue
  if f.resolve()==output:continue
  if f.is_symlink():raise RuntimeError(f'Unexpected backup symlink: {rel}; explicitly scope the backed-up directory')
  if f.is_file():expected[prefix+str(rel)]=f.stat().st_size
 assert expected,'No files selected for backup verification'
 listing=json.loads(subprocess.check_output(aws+['s3api','list-objects-v2','--bucket',BUCKET,'--prefix',prefix,'--query','Contents[].{key:Key,bytes:Size}','--output','json'],text=True)) or []
 actual={x['key']:x['bytes'] for x in listing};assert len(actual)==len(listing)
 missing,mismatched=compare(expected,actual)
 result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'account':account,'bucket':BUCKET,'prefix':prefix,'local_root':str(root),'excluded_directories':a.exclude_dir,'expected_files':len(expected),'expected_bytes':sum(expected.values()),'listed_s3_objects':len(actual),'listed_s3_bytes':sum(actual.values()),'missing':missing,'size_mismatches':mismatched,'verified':not missing and not mismatched,'verification_scope':'Every selected local regular file has an S3 object of equal size. This inventory check is not a cryptographic comparison of object contents. Upload success, source hashes and scientific parity checks are retained separately.'}
 output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2));print(json.dumps(result))
 assert result['verified'],'Private backup inventory mismatch; retain the worker and repair the backup'

if __name__=='__main__':main()

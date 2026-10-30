"""Check final archive symlinks against preserved target entries and exclusions."""
import datetime,hashlib,json,os,subprocess
from pathlib import Path
p=Path(__file__).resolve().parent;closure=json.loads((p/'primary-final-evidence/complete.json').read_text());entries={};exclusions={}
for name in closure['scope_receipts']:
 receipt=json.loads((p/'primary-final-evidence'/name).read_text())
 assert receipt['archive']['independent_download_hash_verified']
 entries.update({x['path']:x for x in receipt['files']})
 exclusions.update({x['path']:x['reason'] for x in receipt['exclusions']})
links=[x for x in entries.values() if 'symlink' in x]
# Enumerate targets without hashing the large datasets a second time. The archive
# manifests already retain the source hashes and independent download verification.
script='''import os,json,pathlib
r=pathlib.Path('/home/ec2-user/experiment');out=[]
for x in json.loads(INPUT):
 target=pathlib.Path(x['resolved_target']);assert target.exists(),x
 if target.is_file():paths=[str(target.relative_to(r))]
 else:
  paths=[]
  for base,dirs,files in os.walk(target,followlinks=False):
   for name in list(dirs):
    q=pathlib.Path(base)/name
    if q.is_symlink():dirs.remove(name);paths.append(str(q.relative_to(r)))
   paths.extend(str((pathlib.Path(base)/name).relative_to(r)) for name in files)
 out.append({'path':x['path'],'target':str(target.relative_to(r)),'target_files':paths})
print(json.dumps(out))
'''.replace('INPUT',repr(json.dumps(links)))
ssh=p.parent/'aa-experiment-20261021/ssh_config'
rows=json.loads(subprocess.check_output(['ssh','-F',str(ssh),'aa-experiment',"python3 - <<'P'\n"+script+"\nP"],text=True));missing=[];covered=[]
for row in rows:
 included=[];excluded=[]
 for target in row['target_files']:
  if target in entries:included.append(target);continue
  reason=next((reason for base,reason in exclusions.items() if target==base or target.startswith(base+'/')),None)
  if reason:excluded.append({'path':target,'reason':reason})
  else:missing.append({'link':row['path'],'target':target})
 covered.append({'path':row['path'],'target':row['target'],'included_target_files':len(included),'documented_excluded_target_files':excluded})
assert not missing,missing[:10]
result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'closure_sha256':hashlib.sha256((p/'primary-final-evidence/complete.json').read_bytes()).hexdigest(),'links':len(links),'coverage':covered,'unpreserved_existing_scientific_targets':missing,'limits':['Transient work and cache exclusions retain their explicit archive reasons; no scientific target is implicitly excused.','File contents are verified by the source-hashed independently downloaded archives, not rehashed by this path-coverage check.']}
(p/'primary-symlink-coverage.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'verified_links':len(links),'unpreserved_existing_scientific_targets':len(missing),'documented_excluded_target_files':sum(len(x['documented_excluded_target_files']) for x in covered)}))

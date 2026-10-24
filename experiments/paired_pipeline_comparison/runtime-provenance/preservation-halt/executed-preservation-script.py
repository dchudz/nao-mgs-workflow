import concurrent.futures,json,pathlib,subprocess
L=pathlib.Path('/Users/david/Documents/ChatGPT/pipeline/.tmp/paired-comparison-20261024')
script=r'''
import datetime,hashlib,json,pathlib,subprocess,sys,tarfile
worker=sys.argv[1];root=pathlib.Path('/home/ec2-user/paired-study-20261024')
aws=['aws','--profile','pipeline-codex','--region','us-east-1']
assert json.loads(subprocess.check_output(aws+['sts','get-caller-identity'],text=True))['Account']=='839993661825'
out=root/'preservation-halt';out.mkdir(exist_ok=True)
units=[]
for d in sorted((root/'comparisons').iterdir()):
 if d.is_dir():
  proof=d/'completed-evidence-backup.json'
  units.append(dict(label=d.name,complete=(d/'complete.json').exists(),backup=json.loads(proof.read_text()) if proof.exists() else None))
inventory=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),worker=worker,units=units,unfinished_state='Retained on EBS; processes suspended; no new analysis requested',mounts=subprocess.check_output(['findmnt','-J'],text=True))
(out/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
archive=out/(worker+'-halt-metadata.tar.gz')
paths=set((root/'batch-logs').rglob('*'))|set((root/'comparisons').rglob('*.json'))|set((root/'comparisons').rglob('*.log'))|{out/'inventory.json',root/'execution-bound.json'}
with tarfile.open(archive,'w:gz',compresslevel=1) as t:
 for f in sorted(paths):
  if f.is_file() and not f.is_symlink() and f.stat().st_size<=10*1024**2:t.add(f,arcname=str(f.relative_to(root)))
key='s3://pipeline-lab-839993661825-us-east-1/runs/paired-comparison-20261024/preservation-halt/'+archive.name
subprocess.run(aws+['s3','cp',str(archive),key,'--only-show-errors'],check=True)
expected=hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()
p=subprocess.Popen(aws+['s3','cp',key,'-','--only-show-errors'],stdout=subprocess.PIPE)
h=hashlib.sha256();size=0
for b in iter(lambda:p.stdout.read(8*1024**2),b''):h.update(b);size+=len(b)
assert p.wait()==0 and size==archive.stat().st_size and h.hexdigest()==expected
receipt=dict(s3=key,bytes=size,sha256=expected,readback_verified=True,inventory=inventory)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
subprocess.run(['sync'],check=True)
print(json.dumps(receipt))
'''
def run(pair):
 host,worker=pair
 r=subprocess.run(['ssh','-F',str(L/'ssh_config'),host,'python3.12 - '+worker],input=script,text=True,capture_output=True,check=True)
 d=json.loads(r.stdout);(L/(worker+'-halt-metadata.verified.json')).write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps({'worker':worker,'verified':True,'bytes':d['bytes'],'units':len(d['inventory']['units']),'preserved_complete':sum(bool(x['backup']) for x in d['inventory']['units'])}),flush=True)
with concurrent.futures.ThreadPoolExecutor(2) as pool:list(pool.map(run,[('paired-study','primary'),('paired-spot2','spot2')]))

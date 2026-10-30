"""Preserve interrupted diagnostics, regenerate only uncommitted generated scratch."""
import datetime,json,shutil,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');prefix='runs/amino-acid-20261021/cohort-screens/'
keys=json.loads(subprocess.check_output(['aws','--profile','pipeline-codex','s3api','list-objects-v2','--bucket','pipeline-lab-839993661825-us-east-1','--prefix',prefix,'--query','Contents[?ends_with(Key, `/complete.json`)].Key','--output','json'],text=True))
committed={x.split('/')[-2] for x in keys};rows=json.loads((r/'heldout-screen-first800.json').read_text());remaining=[x for x in rows if x['run_accession'] not in committed]
records=[];archive=r/'screen-resize-interruption';archive.mkdir(exist_ok=True)
for row in remaining:
 d=r/'cohort/heldout'/row['run_accession']
 if not d.exists() or (d/'complete.json').exists():continue
 saved=archive/d.name;saved.mkdir(exist_ok=True)
 for p in d.iterdir():
  if p.is_file() and p.suffix in ['.json','.log']:shutil.copy2(p,saved/p.name)
 removed={}
 for name in ['reads.fastq','nt.fastq','nt20.fastq','nt22.fastq','aa12.fastq','union.fastq','rescued.fastq']:
  p=d/name
  if p.exists():removed[name]=p.stat().st_size;p.unlink()
 for p in d.glob('*.partial'):removed[p.name]=p.stat().st_size;p.unlink()
 if (d/'failure.json').exists():(d/'failure.json').unlink()
 records.append({'sample':d.name,'regenerated_scratch_bytes':removed})
(archive/'record.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'committed_first800':800-len(remaining),'remaining_first800':len(remaining),'interrupted':records},indent=2))
(r/'heldout-screen-remaining.json').write_text(json.dumps(remaining,indent=2))
subprocess.run(['aws','--profile','pipeline-codex','s3','cp',str(archive)+'/', 's3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/screen-resize-interruption/','--recursive','--only-show-errors'],check=True)
print(json.dumps({'already_committed':800-len(remaining),'remaining':len(remaining),'regenerated':len(records)}))

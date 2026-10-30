"""Compact read-only progress and worker-time bounds for the private experiment."""
import csv,datetime,json,shutil,subprocess
from pathlib import Path
R=Path('/home/ec2-user/experiment');bound=json.loads((R/'execution-bound.json').read_text())
now=datetime.datetime.now(datetime.timezone.utc);hours=(now-datetime.datetime.fromisoformat(bound['prepared_utc'])).total_seconds()/3600
samples={}
for d in sorted((R/'e2e').glob('*')):
    if not d.is_dir() or d.name=='analysis':continue
    complete=(d/'complete.json').exists()
    stages={}
    for trace in sorted(d.glob('*/trace.tsv')):
        rows=list(csv.DictReader(trace.open(),delimiter='\t'))
        stages[trace.parent.name]={'finished_processes':len(rows),'failed_processes':sum(r.get('status')=='FAILED' for r in rows)}
    samples[d.name]={'complete':complete,'stages':stages}
cache=R/'scratch/blast_db_190b4ae1';cached=sum(p.stat().st_size for p in cache.glob('*') if p.is_file())
print(json.dumps({'utc':now.isoformat(),'worker_hours':round(hours,3),'worker_cost_upper_USD':round(hours*4,2),
 'free_GiB':round(shutil.disk_usage(R).free/1024**3,1),'core_cache_GiB':round(cached/1024**3,1),
 'protein_fixture_verified':(R/'integration/protein-workflow/verification.json').exists(),
 'adapter_parity_verified':(R/'integration/adapter-bounded/verification.json').exists(),'samples':samples}))

"""Collect completed search durations without substituting client CPU or incomplete logs."""
import json,re,shlex
from pathlib import Path
r=Path(__file__).resolve().parent/'results';rows=[]
primary_hour=3.3345+(2048*.08+3000*.005+875*.04)/730
nr_hour=2.47752+(1280*.08+3000*.005+875*.04)/730
for task in ['megablast','blastn']:
 p=r/'seed-diagnostic200'/f'{task}.resources.json'
 if p.exists():
  x=json.loads(p.read_text());rows.append({'method':task,'scope':'Same-read diagnostic,746pairs','queries':1436,'wall_seconds':x['wall_seconds'],'cpu_seconds':None,'allocated_worker_cost_usd':x['wall_seconds']/3600*primary_hour,'command':shlex.join(x['command']),'runtime_scope':'Sharedworker; sensitiveBLAST CPUshares changed during execution. Wholeworker allocation is not incremental billing.'})
for folder,method,scope in [('protein-nr','DIAMOND blastx against full NR, pooled','Pooled candidates, native controls, and same-read diagnostic'),('seed-nr-standalone','DIAMOND blastx against full NR, standalone','Same-read diagnostic,746 pairs')]:
 p=r/folder/'search-protocol.json';log=r/folder/'nr-search.log'
 if p.exists() and log.exists():
  text=log.read_text()
  if re.search(r'Exit status:\s*0\s*$',text) and 'Command terminated by signal' not in text:
   x=json.loads(p.read_text());m=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)',text);assert m
   seconds=0
   for v in m[1].split(':'):seconds=seconds*60+float(v)
   cpu=sum(float(re.search(re.escape(label)+r':\s*([\d.]+)',text)[1]) for label in ['User time (seconds)','System time (seconds)'])
   rows.append({'method':method,'scope':scope,'queries':x['queries'],'wall_seconds':seconds,'cpu_seconds':cpu,'allocated_worker_cost_usd':seconds/3600*nr_hour,'command':shlex.join(x['command']),'runtime_scope':'Completed full-reference scan and per-search seed indexing. Pooled runtime is not allocated to subsets. Hardware and scheduling differ from nucleotide runs.'})
(r/'search-runtimes.json').write_text(json.dumps({'runs':rows,'notes':['Unmeasured CPU remainsnull. Wholeworker cost allocations must not be summed across concurrent jobs.','Interrupted attempts, reference staging, and idle time belong in the overall instance ledger.']},indent=2))
print(json.dumps([{k:v for k,v in x.items() if k not in ['command','runtime_scope']} for x in rows]))

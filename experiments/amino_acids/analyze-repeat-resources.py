"""Measure paired downstream CPU from existing traces; do not infer full-raw repeat-screen speed."""
import collections,csv,json,re,statistics
from pathlib import Path
r=Path('/home/ec2-user/experiment');c=r/'cohort-storage';d=c/'repeat-confirmation';out=d/'analysis'
planned=set(json.loads((r/'repeat-confirmation-protocol.json').read_text())['samples']);assert len(planned)==500
original={x['sample']:x for x in csv.DictReader((c/'analysis/resource-metrics.tsv').open(),delimiter='\t') if x['arm']=='union' and x['sample'] in planned};assert original.keys()==planned
units={'ms':.001,'s':1,'m':60,'h':3600,'d':86400}
def seconds(value):
 if value=='-':return None
 parts=re.findall(r'([\d.]+)(ms|s|m|h|d)',value);assert parts,value
 return sum(float(n)*units[u] for n,u in parts)
tasks=collections.defaultdict(list);screens={};completed=[]
for marker in sorted((d/'batches').glob('*/complete.json')):
 batch=marker.parent;completed.extend(json.loads(marker.read_text())['samples'])
 for x in json.loads((batch/'screen-metrics.json').read_text()):
  assert x['sample'] not in screens;screens[x['sample']]=x
 for x in csv.DictReader((batch/'output/trace.tsv').open(),delimiter='\t'):
  m=re.search(r'(SRR\d+)_repeatunion\b',x['name']);assert m,x['name']
  elapsed=seconds(x['realtime']);cpu=None if x['%cpu']=='-' else float(x['%cpu'].rstrip('%'));tasks[m[1]].append(None if elapsed is None or cpu is None else elapsed*cpu/100)
assert len(completed)==len(set(completed))==500 and set(completed)==planned==screens.keys()==tasks.keys()
rows=[]
for sample in sorted(planned):
 old=original[sample];new=tasks[sample];screen=screens[sample]['subset_screen_resources']
 row={'sample':sample,'original_union_downstream_CPU_seconds_estimate':float(old['downstream_CPU_seconds_estimate']),'filtered_union_downstream_CPU_seconds_estimate':sum(x for x in new if x is not None),'original_tasks_missing_CPU':int(old['trace_tasks_missing_CPU']),'filtered_tasks_missing_CPU':sum(x is None for x in new),'filtered_tasks':len(new),'subset_rescreen_CPU_seconds':screen['user_seconds']+screen['system_seconds']}
 row['downstream_CPU_change_seconds_estimate']=row['filtered_union_downstream_CPU_seconds_estimate']-row['original_union_downstream_CPU_seconds_estimate'];rows.append(row)
with (out/'resource-metrics.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
metrics={k:{'total':sum(x[k] for x in rows),'median_per_sample':statistics.median(x[k] for x in rows)} for k in rows[0] if k!='sample'}
summary={'samples':500,'metrics':metrics,'all_tasks_have_CPU':all(x['original_tasks_missing_CPU']==x['filtered_tasks_missing_CPU']==0 for x in rows),'aggregate_downstream_CPU_reduction_fraction':1-metrics['filtered_union_downstream_CPU_seconds_estimate']['total']/metrics['original_union_downstream_CPU_seconds_estimate']['total'],'interpretation':['Both downstream measurements use Nextflow realtime multiplied by task percent CPU; they are estimates under differing concurrent loads.','The correction was screened over already retained AA12 candidates, not full raw inputs. Its GNU-time CPU is separately reported and cannot establish full-pipeline speed.','Excludes extraction, index construction, union construction, orchestration, storage and idle workers.']}
(out/'resource-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

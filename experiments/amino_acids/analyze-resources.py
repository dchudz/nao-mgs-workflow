"""Stage-wise resource comparison; distinguish measured screen CPU from trace estimates."""
import collections,csv,json,re,statistics,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/analysis';batch_root=Path(sys.argv[1]) if len(sys.argv)>1 else r/'cohort-storage/batches';arms=['nt','nt22','nt20','aa12','union']
def seconds(value):
 if value=='-':return None
 units={'ms':.001,'s':1,'m':60,'h':3600,'d':86400};parts=re.findall(r'([\d.]+)(ms|s|m|h|d)',value)
 assert parts,value
 return sum(float(n)*units[u] for n,u in parts)
counts={'nt':'nucleotide_pairs','aa12':'amino_pairs','union':'union_pairs'}
screen={};processed=collections.defaultdict(lambda:collections.defaultdict(list))
for marker in batch_root.glob('*/complete.json'):
 b=marker.parent
 for a in json.loads(marker.read_text())['samples']:
  p=r/'cohort-storage/screen-burden-inputs'/a/'complete.json'
  assert p.exists(), 'Missing measured screen resources: '+a
  screen[a]=json.loads(p.read_text())
 with open(b/'output/trace.tsv') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   m=re.search(r'(SRR\d+)_(nt22|nt20|aa12|union|nt)\b',x['name'])
   if not m:raise ValueError('Unattributed task: '+x['name'])
   a,arm=m.groups();elapsed=seconds(x['realtime']);cpu=None if x['%cpu']=='-' else float(x['%cpu'].rstrip('%'))
   processed[a][arm].append(None if elapsed is None or cpu is None else elapsed*cpu/100)
rows=[]
for a,arms_data in processed.items():
 assert a in screen, 'Missing measured screen resources: '+a
 s=screen[a]
 for arm in arms:
  resources=s['resources'];selected=['nt','aa12'] if arm=='union' else [arm];screen_cpu=sum(resources[k]['user_seconds']+resources[k]['system_seconds'] for k in selected)
  tasks=arms_data[arm];row={'sample':a,'arm':arm,'screen_CPU_seconds':screen_cpu,'downstream_CPU_seconds_estimate':sum(x for x in tasks if x is not None),'trace_tasks':len(tasks),'trace_tasks_missing_CPU':sum(x is None for x in tasks),'screen_peak_rss_kib':max(resources[k]['peak_rss_kib'] for k in selected)}
  row['screen_selected_pairs']=s.get(counts.get(arm,''))
  rows.append(row)
with (d/'resource-metrics.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['sample'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'samples':len(rows)//5,'arms':{},'interpretation':'Screen CPU is native GNU time user+system. Downstream CPU is an estimate from Nextflow realtime times percent-CPU. Excludes extraction, index construction, union construction, orchestration, reference staging, and idle worker time. Wall times under shared concurrency are not isolated speed benchmarks. Downstream hardware and concurrency differ between batches: first33samples and recovery batches use the primary worker; intervening75samples use Spot. Each sample compares all5arms on the same worker.'}
for arm in arms:
 selected=[x for x in rows if x['arm']==arm]
 if selected:summary['arms'][arm]={k:{'total':sum(x[k] for x in selected),'median_per_sample':statistics.median(x[k] for x in selected)} for k in ['screen_CPU_seconds','downstream_CPU_seconds_estimate','trace_tasks_missing_CPU']}
(d/'resource-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
# Full-depth runs are separate: NT24 CPU measured only the Docker client and NT22 screened an NT20 subset.
p=r/'full/comparison/output/trace.tsv';full=collections.defaultdict(float)
for x in csv.DictReader(p.open(),delimiter='\t'):
 m=re.search(r'SRR\d+_(nt22|nt20|aa12|union|nt)\b',x['name']);assert m
 elapsed=seconds(x['realtime']);cpu=float(x['%cpu'].rstrip('%'));full[m[1]]+=elapsed*cpu/100
(d/'full-depth-downstream-CPU.json').write_text(json.dumps({'estimated_CPU_seconds':dict(full),'interpretation':'Downstream trace estimate only. Full-depthNT24 screen CPU is unavailable because GNU time measured Docker client. Full-depthNT22 timing is not comparable because it screened theNT20 subset.'},indent=2))

"""Describe sample-species assignment changes using original pipeline clade counts."""
import collections,csv,gzip,json,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');out=r/'cohort-storage/analysis';batch_root=Path(sys.argv[1]) if len(sys.argv)>1 else r/'cohort-storage/batches';arms=['nt','nt22','nt20','aa12','union']
with gzip.open(r/'reference/results/total-virus-db-annotated.tsv.gz','rt') as f:species={x['taxid']:x for x in csv.DictReader(f,delimiter='\t') if x['rank']=='species'}
changes=[];totals=collections.defaultdict(collections.Counter);samples=0
for marker in batch_root.glob('*/complete.json'):
 b=marker.parent
 for a in json.loads(marker.read_text())['samples']:
  samples+=1;counts={}
  for arm in arms:
   name=f'{a}_{arm}_clade_counts.tsv.gz';p=b/'output/counts'/name
   if not p.exists():p=r/'cohort-storage/counts-backfill'/name
   with gzip.open(p,'rt') as f:counts[arm]={x['taxid']:int(x['reads_clade_dedup']) for x in csv.DictReader(f,delimiter='\t') if x['taxid'] in species}
  for minimum in [1,2,3]:
   base={t for t,n in counts['nt'].items() if n>=minimum}
   for arm in arms:
    called={t for t,n in counts[arm].items() if n>=minimum};gained=called-base;lost=base-called;key=f'{arm}/minimum_{minimum}'
    totals[key].update({'sample_species_calls':len(called),'gained_sample_species_calls':len(gained),'lost_sample_species_calls':len(lost),'human_species_calls':sum(species[t]['infection_status_human']=='1' for t in called),'gained_human_species_calls':sum(species[t]['infection_status_human']=='1' for t in gained),'lost_human_species_calls':sum(species[t]['infection_status_human']=='1' for t in lost)})
    for t in sorted(gained|lost):changes.append({'sample':a,'arm':arm,'minimum_exemplars':minimum,'change':'gained' if t in gained else 'lost','taxid':t,'name':species[t]['name'],'human_annotation':species[t]['infection_status_human'],'baseline_exemplars':counts['nt'].get(t,0),'arm_exemplars':counts[arm].get(t,0)})
with (out/'species-changes.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(changes[0]) if changes else ['sample'],delimiter='\t');w.writeheader();w.writerows(changes)
summary={'samples':samples,'thresholds':[1,2,3],'counts':{k:dict(v) for k,v in totals.items()},'interpretation':'Original pipeline reads_clade_dedup at species-rank nodes. Thresholds are descriptive support levels, not validated diagnostic cutoffs. Species labels and human-infection annotations come from the fixed reference. A gained sample-species assignment does not establish a new infection or confirmed biological detection.'};(out/'species-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

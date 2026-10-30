"""Secondary paired uncertainty clustered by public BioSample location labels."""
import collections,csv,json,random,statistics,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/analysis'
meta={x['run_accession']:x for x in csv.DictReader(open(r/'heldout-biosample-metadata.tsv'),delimiter='\t')}
rows=list(csv.DictReader(open(d/'sample-metrics.tsv'),delimiter='\t'));result={};site_rows=[];month_rows=[]
for arm in ['nt22','nt20','aa12','union']:
 selected=[x for x in rows if x['arm']==arm];groups=collections.defaultdict(list);months=collections.defaultdict(list)
 for x in selected:
  m=meta[x['sample']];assert m['geo_loc_name'];gain=int(x['added_pairs'])-int(x['lost_pairs']);groups[m['geo_loc_name']].append(gain);months[m['collection_date'][:7]].append(gain)
 rng=random.Random(20261021);keys=sorted(groups);means=[]
 if not keys:continue
 for _ in range(5000):
  sample=rng.choices(keys,k=len(keys));means.append(sum(sum(groups[g]) for g in sample)/sum(len(groups[g]) for g in sample))
 means.sort();result[arm]={'samples':len(selected),'location_clusters':len(keys),'mean_net_gain_per_sample':sum(sum(x) for x in groups.values())/len(selected),'location_cluster_bootstrap_95':[means[125],means[4874]]}
 for g,values in sorted(groups.items()):site_rows.append({'arm':arm,'location_label':g,'samples':len(values),'net_added_pairs':sum(values),'mean_net_gain':statistics.mean(values)})
 for month,values in sorted(months.items()):month_rows.append({'arm':arm,'month':month,'samples':len(values),'net_added_pairs':sum(values),'mean_net_gain':statistics.mean(values)})
for name,data in [('site-effects.tsv',site_rows),('month-effects.tsv',month_rows)]:
 with (d/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['arm'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'secondary':True,'selection_unchanged':True,'arms':result,'interpretation':'Resample whole location-label clusters, retaining their repeated samples and times. Ratio estimator weights the observed samples within selected clusters. Intervals are descriptive for this cohort; location labels can share a city, and the cohort is not a representative population survey.'};(d/'clustered-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

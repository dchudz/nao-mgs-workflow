"""Compare aggregate support counts within the prespecified nested-depth pairs."""
import argparse,csv,hashlib,json,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--audit',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
cases={p.parent.name:json.loads(p.read_text()) for p in a.audit.glob('*/complete.json')}
paired=sorted(k[:-4] for k in cases if k.endswith('_25m') and k[:-4]+'_5m' in cases)
assert len(paired)==10
rows=[];transitions=[];proof={}
for sample in paired:
 for arm in ('production','protein_variant'):
  lo=cases[sample+'_5m']['arms'][arm];hi=cases[sample+'_25m']['arms'][arm]
  for n in (2,3,5):
   low={k for k,v in lo['species_exemplars'].items() if v>=n};high={k for k,v in hi['species_exemplars'].items() if v>=n}
   rows.append(dict(sample=sample,arm=arm,threshold=n,shallow_calls=len(low),deep_calls=len(high),retained=len(low&high),deep_only=len(high-low),shallow_only=len(low-high),net_calls=len(high)-len(low)))
 for depth in (5,25):
  path=a.audit/f'{sample}_{depth}m'/'complete.json';proof[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
summary=[]
for arm in ('production','protein_variant'):
 for n in (2,3,5):
  subset=[r for r in rows if r['arm']==arm and r['threshold']==n]
  r=dict(arm=arm,threshold=n,paired_runs=len(subset))
  for k in ('shallow_calls','deep_calls','retained','deep_only','shallow_only','net_calls'):r[k]=sum(x[k] for x in subset)
  r['median_net_per_run']=statistics.median(x['net_calls'] for x in subset)
  for direction in ('positive','zero','negative'):r['runs_with_'+direction+'_net']=sum((x['net_calls']>0 if direction=='positive' else x['net_calls']==0 if direction=='zero' else x['net_calls']<0) for x in subset)
  summary.append(r)
for name,data in [('depth-per-run',rows),('depth-summary',summary)]:
 with (a.out/(name+'.tsv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
(a.out/'depth-provenance.json').write_text(json.dumps(dict(input_hashes=proof,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Counts from existing outputs only; no new alignments, taxonomic reassignment, or validation searches.'),indent=2)+'\n')
print(json.dumps(summary))

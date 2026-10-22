"""Count paired screen candidates from frozen native logs and AA union records."""
import csv,json,re,statistics,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=Path(sys.argv[1]) if len(sys.argv)>1 else r/'cohort/heldout';out=Path(sys.argv[2]) if len(sys.argv)>2 else r/'screen-burden';out.mkdir(exist_ok=True);rows=[]
for marker in sorted(d.glob('*/complete.json')):
 x=json.loads(marker.read_text());sample=x['run_accession'];counts={}
 for arm in ['nt','nt22','nt20']:
  m=re.search(r'^Matches:\s+(\d+) reads', (marker.parent/(arm+'.log')).read_text(),re.M);assert m,(sample,arm)
  reads=int(m[1]);assert reads%2==0;counts[arm]=reads//2
 assert counts['nt']==x['nucleotide_pairs']
 counts.update(aa12=x['amino_pairs'],union=x['union_pairs'])
 for arm,count in counts.items():rows.append(dict(sample=sample,arm=arm,input_pairs=x['observed_pairs'],candidate_pairs=count,candidate_fraction=count/x['observed_pairs'],candidate_ratio_to_nt24=count/counts['nt'] if counts['nt'] else None))
with (out/'sample-metrics.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={'samples':len(rows)//5,'arms':{},'scope':'Only completed raw screens; all arms see identical raw pairs per sample. Counts precede FASTP and downstream filtering.'}
for arm in ['nt','nt22','nt20','aa12','union']:
 s=[x for x in rows if x['arm']==arm];summary['arms'][arm]={'total_candidate_pairs':sum(x['candidate_pairs'] for x in s),'median_candidate_pairs':statistics.median(x['candidate_pairs'] for x in s),'median_candidate_fraction':statistics.median(x['candidate_fraction'] for x in s),'median_within_sample_ratio_to_nt24':statistics.median(x['candidate_ratio_to_nt24'] for x in s if x['candidate_ratio_to_nt24'] is not None)}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

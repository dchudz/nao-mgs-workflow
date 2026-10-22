"""Paired, sample-level comparisons; only successfully completed batches enter."""
import argparse,collections,csv,gzip,json,random,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('output');a=p.parse_args()
root=Path(a.root);out=Path(a.output);out.mkdir(exist_ok=True)
arms=['nt','nt22','nt20','aa12','union'];rows=[];gains=[]
def read(path):
 with gzip.open(path,'rt') as f:
  rows=list(csv.DictReader(f,delimiter='\t'))
  result={x['seq_id']:x for x in rows};assert len(rows)==len(result),path
  return result
for marker in sorted(root.glob('*/complete.json')):
 batch=marker.parent
 for sample in json.loads(marker.read_text())['samples']:
  data={}
  for arm in arms:
   files=list((batch/'output/marked').glob(f'*_{sample}_{arm}_duplicate_reads.tsv.gz'));assert len(files)==1,(sample,arm)
   data[arm]=read(files[0])
  baseline=data['nt'];baseids=set(baseline)
  for arm in arms:
   d=data[arm];ids=set(d);extra=ids-baseids;lost=baseids-ids
   row=dict(sample=sample,arm=arm,final_pairs=len(d),added_pairs=len(extra),lost_pairs=len(lost),stored_sequence_signatures=len({(x['query_seq'],x['query_seq_rev']) for x in d.values()}),joint_exemplars=sum(x['seq_id']==x['prim_align_dup_exemplar']==x['sim_dup_exemplar'] for x in d.values()),taxid_changes=sum(d[q]['aligner_taxid_lca']!=baseline[q]['aligner_taxid_lca'] for q in ids&baseids),added_beyond_nt20=len(extra-set(data['nt20'])),added_beyond_aa12_union=len(extra-set(data['union'])))
   rows.append(row)
   for q in sorted(extra):gains.append(dict(sample=sample,arm=arm,seq_id=q,taxid=d[q]['aligner_taxid_lca'],joint_exemplar=q==d[q]['prim_align_dup_exemplar']==d[q]['sim_dup_exemplar'],also_nt20=q in data['nt20'],also_aa12_union=q in data['union']))
for name,data in [('sample-metrics.tsv',rows),('added-pairs.tsv',gains)]:
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['sample'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'completed_samples':len(rows)//len(arms),'sampled_input_pairs':len(rows)//len(arms)*5000000,'arms':{},'uncertainty':'Bootstrap intervals resample biological samples, not reads. Related sites/dates remain correlated; intervals are descriptive and may understate uncertainty.'}
rng=random.Random(20261021)
for arm in arms:
 selected=[x for x in rows if x['arm']==arm]
 if not selected:continue
 totals={k:sum(x[k] for x in selected) for k in selected[0] if k not in ['sample','arm']}
 deltas=[x['added_pairs']-x['lost_pairs'] for x in selected]
 bootstrap=sorted(sum(rng.choices(deltas,k=len(deltas)))/len(deltas) for _ in range(2000))
 summary['arms'][arm]=dict(totals=totals,samples_with_gain=sum(x['added_pairs']>0 for x in selected),samples_with_loss=sum(x['lost_pairs']>0 for x in selected),median_net_gain=statistics.median(deltas),mean_net_gain=sum(deltas)/len(deltas),mean_gain_bootstrap_95=[bootstrap[50],bootstrap[1949]])
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

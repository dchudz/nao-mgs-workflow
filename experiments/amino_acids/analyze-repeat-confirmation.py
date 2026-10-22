"""Paired prospective last-500 comparison; report every gain and loss, including native losses."""
import collections,csv,gzip,json,statistics
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/repeat-confirmation';out=d/'analysis';out.mkdir(exist_ok=True)
original={s:m.parent for m in (r/'cohort-storage/batches').glob('*/complete.json') for s in json.loads(m.read_text())['samples']}
planned=json.loads((r/'repeat-confirmation-protocol.json').read_text())['samples'];assert len(planned)==500
rows=[];changes=[];pending=[]
def read(p):
 with gzip.open(p,'rt') as f:
  a=list(csv.DictReader(f,delimiter='\t'));z={x['seq_id']:x for x in a};assert len(a)==len(z);return z
for marker in (d/'batches').glob('*/complete.json'):
 b=marker.parent;screens={x['sample']:x for x in json.loads((b/'screen-metrics.json').read_text())}
 for s in json.loads(marker.read_text())['samples']:
  assert s in planned
  if s not in original:pending.append(s);continue
  tables={}
  for arm,base in [('nt',original[s]),('union',original[s]),('repeatunion',b)]:
   paths=list((base/'output/marked').glob(f'*_{s}_{arm}_duplicate_reads.tsv.gz'));assert len(paths)==1
   tables[arm]=read(paths[0])
  native,old,new=(tables[a] for a in ['nt','union','repeatunion']);extra=new.keys()-old.keys();lost=old.keys()-new.keys();v=screens[s]
  row=dict(sample=s,original_candidate_pairs=v['original_union_pairs'],repeat_candidate_pairs=v['union_pairs'],candidate_reduction_fraction=1-v['union_pairs']/v['original_union_pairs'],native_final_pairs=len(native),original_union_final_pairs=len(old),repeat_final_pairs=len(new),lost_original_union_pairs=len(lost),new_vs_original_union_pairs=len(extra),lost_native_pairs=len(native.keys()-new.keys()),added_to_native_pairs=len(new.keys()-native.keys()),taxid_changes=sum(new[q]['aligner_taxid_lca']!=old[q]['aligner_taxid_lca'] for q in new.keys()&old.keys()),original_joint_exemplars=sum(x['seq_id']==x['prim_align_dup_exemplar']==x['sim_dup_exemplar'] for x in old.values()),repeat_joint_exemplars=sum(x['seq_id']==x['prim_align_dup_exemplar']==x['sim_dup_exemplar'] for x in new.values()))
  common=old.keys()&new.keys()
  row['shared_query_sequence_changes']=sum((old[q]['query_seq'],old[q]['query_seq_rev'])!=(new[q]['query_seq'],new[q]['query_seq_rev']) for q in common)
  row['primary_duplicate_exemplar_changes']=sum(old[q]['prim_align_dup_exemplar']!=new[q]['prim_align_dup_exemplar'] for q in common)
  row['similarity_duplicate_exemplar_changes']=sum(old[q]['sim_dup_exemplar']!=new[q]['sim_dup_exemplar'] for q in common)
  row['joint_exemplar_membership_changes']=sum((old[q]['seq_id']==old[q]['prim_align_dup_exemplar']==old[q]['sim_dup_exemplar'])!=(new[q]['seq_id']==new[q]['prim_align_dup_exemplar']==new[q]['sim_dup_exemplar']) for q in common)
  rows.append(row)
  for q in sorted(extra|lost):changes.append(dict(sample=s,seq_id=q,change='gained' if q in extra else 'lost',taxid=(new if q in extra else old)[q]['aligner_taxid_lca'],also_native=q in native))
for name,data in [('sample-metrics.tsv',rows),('changed-pairs.tsv',changes)]:
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['sample'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'secondary_prospective':True,'completed_paired_samples':len(rows),'planned_samples':500,'repeat_completed_waiting_for_primary':pending,'rule':'Omit AA12 words with ten or more of the same residue; preserve original NT24 candidates.','totals':{k:sum(x[k] for x in rows) for k in rows[0] if k not in ['sample','candidate_reduction_fraction']} if rows else {},'median_candidate_reduction_fraction':statistics.median(x['candidate_reduction_fraction'] for x in rows) if rows else None,'samples_with_original_union_loss':sum(x['lost_original_union_pairs']>0 for x in rows),'samples_with_native_loss':sum(x['lost_native_pairs']>0 for x in rows),'interpretation':'Rule selected from an earlier outlier and frozen before any of these500samples started. This is a secondary prospective comparison; the primary1000-sample parameters remain unchanged. Assignment preservation is not known-truth biological sensitivity.'}
if rows:summary['aggregate_candidate_reduction_fraction']=1-summary['totals']['repeat_candidate_pairs']/summary['totals']['original_candidate_pairs']
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

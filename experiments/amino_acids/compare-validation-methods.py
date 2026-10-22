"""Join three methods on the same original read pairs; retain every disagreement."""
import argparse,collections,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--nucleotide',required=True);p.add_argument('--protein',required=True);p.add_argument('--output',required=True);p.add_argument('--runtime-scope',choices=['pooled','standalone'],default='pooled');a=p.parse_args()
nt=list(csv.DictReader(open(a.nucleotide),delimiter='\t'));protein=list(csv.DictReader(open(a.protein),delimiter='\t'))
by={}
for x in nt:
 k=(x['task'],x['seq_id']);assert k not in by;by[k]=x
nr={x['query_pair'].removeprefix('SEED_'):x for x in protein};assert len(nr)==len(protein)
assert {task for task,q in by}=={'megablast','blastn'}
selection=list(csv.DictReader(open(Path(a.nucleotide).parent/'selection.tsv'),delimiter='\t'))
selected_ids={x['seq_id'] for x in selection};assert len(selected_ids)==len(selection)
ids={q for task,q in by if task=='megablast'}
assert ids==selected_ids,'Comparison omits a selected original read pair'
assert ids=={q for task,q in by if task=='blastn'}==set(nr)
rows=[]
for q in sorted(ids):
 m=by['megablast',q];b=by['blastn',q];n=nr[q]
 assert all(m[k]==b[k] for k in ['seq_id','sample','AA12_added','NT20_added','baseline_control','original_taxid']),'Nucleotide metadata mismatch'
 row={k:m[k] for k in ['seq_id','sample','AA12_added','NT20_added','baseline_control','original_taxid']}
 for name,x,fields in [('megablast',m,['status','lca','aligned_mates']),('blastn',b,['status','lca','aligned_mates']),('protein',n,['nr_status','nr_taxid_lca','aligned_mates'])]:
  for dest,source in zip(['status','lca','aligned_mates'],fields):row[name+'_'+dest]=x[source]
 row['status_disagreement']=len({row[k+'_status'] for k in ['megablast','blastn','protein']})>1
 row['taxid_disagreement']=len({row[k+'_lca'] for k in ['megablast','blastn','protein']})>1
 rows.append(row)
out=Path(a.output);out.mkdir(exist_ok=True,parents=True)
for name,subset in [('all-pairs.tsv',rows),('disagreements.tsv',[x for x in rows if x['status_disagreement'] or x['taxid_disagreement']])]:
 with open(out/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(subset)
summary={'protein_runtime_scope':a.runtime_scope,'pairs':len(rows),'available_query_sequences':sum(int(x['available_mates']) for x in selection),'groups':{},'limitations':['Identical available mate sequences, searched separately; original production validation joins or merges mates.','Nucleotide methods search core_nt; protein method searches NR. Reference coverage and representation differ, so this is not an isolated aligner comparison.',('Protein runtime belongs to the pooled comprehensive-reference job; per-diagnostic runtime is not measured.' if a.runtime_scope=='pooled' else 'Protein runtime is measured by a separate full-NR search over these1,436availablematequeries; hardware and scheduling differ from nucleotide runs.'),'Near-best taxonomy filters differ between production nucleotide ranking and the frozen protein protocol. Disagreement is evidence to inspect, not a known-truth accuracy verdict.']}
for group in ['all','AA12_added','NT20_added','baseline_control']:
 selected=[x for x in rows if group=='all' or x[group]=='True'];z={'pairs':len(selected),'status_disagreements':sum(x['status_disagreement'] for x in selected),'taxid_disagreements':sum(x['taxid_disagreement'] for x in selected),'methods':{}}
 for method in ['megablast','blastn','protein']:z['methods'][method]=dict(collections.Counter(x[method+'_status'] for x in selected))
 for method in ['blastn','protein']:
  z[method+'_viral_only_where_megablast_unaligned']=sum(x['megablast_status']=='no_alignment' and x[method+'_status']=='viral_only' for x in selected)
 summary['groups'][group]=z
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

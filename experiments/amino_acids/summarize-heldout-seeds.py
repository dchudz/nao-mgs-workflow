"""Compare identical separate-mate queries under two BLAST task configurations."""
import collections,csv,json,sys
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=Path(sys.argv[1])
parents={}
for line in open(r/'reference/results/taxonomy-nodes.dmp'):
 a=line.split('|');parents[int(a[0])]=int(a[1])
cache={}
def path(t):
 if t in cache:return cache[t]
 seen=set();result=[];o=t
 while t in parents and t not in seen:
  seen.add(t);result.append(t)
  if t==parents[t]:cache[o]=result;return result
  t=parents[t]
 cache[o]=None;return None
def lca(taxa):
 ps=[path(t) for t in taxa]
 if not ps or any(p is None for p in ps):return 0
 common=set(ps[0]).intersection(*(set(p) for p in ps[1:]));return next((t for t in ps[0] if t in common),0)
selection={x['seq_id']:x for x in csv.DictReader(open(d/'selection.tsv'),delimiter='\t')};out=[]
for task in ['megablast','blastn']:
 hits=collections.defaultdict(list)
 for row in csv.reader(open(d/(task+'.tsv')),delimiter='\t'):
  q,subject,taxa,score,*rest=row;hits[q].append((float(score),{int(t) for t in taxa.split(';') if t.isdigit()} or {0}))
 calls={}
 for q,alignments in hits.items():
  scores=sorted({s for s,t in alignments},reverse=True);threshold=max(0.9*scores[0],scores[min(9,len(scores)-1)])
  calls[q]=set().union(*(t for s,t in alignments if s>=threshold))
 for q,sel in selection.items():
  taxa=set().union(*(calls.get(q+'/'+str(i),set()) for i in [1,2]));unknown={t for t in taxa if path(t) is None};known=taxa-unknown
  viral={t for t in known if 10239 in path(t)};nonviral=known-viral
  artificial={t for t in nonviral if 81077 in path(t)};cellular={t for t in nonviral-artificial if any(root in path(t) for root in [131567,2,2157,2759])};other=nonviral-artificial-cellular
  status='no_alignment' if not taxa else 'mixed_viral_nonviral' if viral and nonviral else 'unknown_taxonomy' if unknown else 'nonviral_only' if nonviral else 'viral_only'
  taxid=lca(taxa);original=int(sel['original_taxid']);a=path(original);b=path(taxid)
  relation='unknown' if a is None or b is None else 'exact' if original==taxid else 'ancestor_related' if original in b or taxid in a else 'divergent'
  out.append(dict(seq_id=q,sample=sel['sample'],AA12_added=sel['AA12_added'],NT20_added=sel['NT20_added'],baseline_control=sel['baseline_control'],task=task,original_taxid=original,lca=taxid,status=status,relation=relation,aligned_mates=sum(q+'/'+str(i) in calls for i in [1,2]),viral_taxa=len(viral),artificial_taxa=len(artificial),cellular_taxa=len(cellular),other_known_taxa=len(other),unknown_taxa=len(unknown),viral_and_artificial_only=bool(viral and artificial and not (cellular or other or unknown))))
with open(d/'comparison.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary={'interpretation':'Heldout first-200 secondary diagnostic with the same available individual mate sequences in both tasks. Original pipeline joins/merges mates, so this is not an exact replay of original validation. Taxonomy uses original top 10 dense scores and 90% near-best filters separately per mate, then unions mate taxa. The task change includes seed and scoring defaults, not an isolated word-size change. Alignment is not known-truth validation.','groups':{}}
for task in ['megablast','blastn']:
 for stratum in ['AA12_added','NT20_added','baseline_control']:
  rows=[x for x in out if x['task']==task and x[stratum]=='True']
  summary['groups'][task+'/'+stratum]={'pairs':len(rows),'samples':len({x['sample'] for x in rows}),'status':dict(collections.Counter(x['status'] for x in rows)),'relation':dict(collections.Counter(x['relation'] for x in rows)),'viral_and_artificial_only':sum(x['viral_and_artificial_only'] for x in rows),'any_cellular_competitor':sum(x['cellular_taxa']>0 for x in rows)}
summary['strata_overlap']=True
summary['selection']='Up to 2 SHA-ranked pairs per sample per stratum; samples with fewer qualifying pairs contribute fewer.'
summary['paired_transitions']={}
by_task={(x['task'],x['seq_id']):x for x in out}
for stratum in ['AA12_added','NT20_added','baseline_control']:
 rows=[x for x in out if x['task']=='megablast' and x[stratum]=='True']
 summary['paired_transitions'][stratum]=dict(collections.Counter(x['status']+' -> '+by_task[('blastn',x['seq_id'])]['status'] for x in rows))
(d/'comparison.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

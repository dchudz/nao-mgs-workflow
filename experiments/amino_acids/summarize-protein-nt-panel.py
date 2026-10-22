"""Keep every panel pair and join observed nucleotide and protein evidence without extrapolation."""
import argparse,collections,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('panel');p.add_argument('--nodes',required=True);p.add_argument('--nr');a=p.parse_args();d=Path(a.panel)
parents={}
for line in open(a.nodes):
 x=line.split('|');parents[int(x[0])]=int(x[1])
cache={}
def ancestry(t):
 if t in cache:return cache[t]
 original=t;seen=set();out=[]
 while t in parents and t not in seen:
  seen.add(t);out.append(t)
  if t==parents[t]:cache[original]=out;return out
  t=parents[t]
 cache[original]=None;return None
def lca(taxa):
 paths=[ancestry(t) for t in taxa]
 if not paths or any(x is None for x in paths):return 0
 common=set(paths[0]).intersection(*(set(x) for x in paths[1:]));return next((t for t in paths[0] if t in common),0)
selection={x['query_pair']:x for x in csv.DictReader((d/'selection.tsv').open(),delimiter='\t')};rows=[];details=[]
for task in ['megablast','blastn']:
 hits=collections.defaultdict(list)
 for x in csv.reader((d/f'{task}.tsv').open(),delimiter='\t'):
  q,subject,taxa,score,evalue,length,identity,coverage=x;assert q.rsplit('/',1)[0] in selection
  hits[q].append(dict(subject=subject,taxids=taxa,bitscore=float(score),evalue=evalue,alignment_nt=int(length),identity_pct=float(identity),coverage_pct=float(coverage)))
 calls={}
 for q,values in hits.items():
  scores=sorted({x['bitscore'] for x in values},reverse=True);minimum=max(.9*scores[0],scores[min(9,len(scores)-1)])
  kept=[x for x in values if x['bitscore']>=minimum]
  calls[q]=set().union(*({int(t) for t in x['taxids'].split(';') if t.isdigit()} or {0} for x in kept))
  for x in kept:details.append(dict(task=task,query=q,**x))
 for q,sel in selection.items():
  taxa=set().union(*(calls.get(q+'/'+str(m),set()) for m in [1,2]));unknown={t for t in taxa if ancestry(t) is None};known=taxa-unknown
  viral={t for t in known if 10239 in ancestry(t)};nonviral=known-viral;artificial={t for t in nonviral if 81077 in ancestry(t)};cellular={t for t in nonviral-artificial if any(x in ancestry(t) for x in [131567,2,2157,2759])}
  status='no_alignment' if not taxa else 'mixed_viral_nonviral' if viral and nonviral else 'unknown_taxonomy' if unknown else 'nonviral_only' if nonviral else 'viral_only'
  taxid=lca(taxa);ref=int(sel['viral_reference_taxid']);pa=ancestry(ref);pb=ancestry(taxid)
  relation='unknown' if pa is None or pb is None else 'exact' if ref==taxid else 'ancestor_related' if ref in pb or taxid in pa else 'divergent'
  rows.append(dict(query_pair=q,sample=sel['sample'],stratum=sel['stratum'],task=task,status=status,taxid_lca=taxid,relation_to_viral_reference=relation,aligned_mates=sum(q+'/'+str(m) in calls for m in [1,2]),viral_taxa=len(viral),artificial_taxa=len(artificial),cellular_taxa=len(cellular),unknown_taxa=len(unknown),near_best_taxids=';'.join(map(str,sorted(taxa)))))
if a.nr:
 nr={x['query_pair']:x for x in csv.DictReader(open(a.nr),delimiter='\t') if x['query_pair'] in selection};assert nr.keys()==selection.keys(),'Missing measured NR evidence'
 for q,sel in selection.items():
  x=nr[q];rows.append(dict(query_pair=q,sample=sel['sample'],stratum=sel['stratum'],task='protein_NR',status=x['nr_status'],taxid_lca=int(x['nr_taxid_lca']),relation_to_viral_reference=x['taxonomic_relation'],aligned_mates=int(x['aligned_mates']),viral_taxa=x['viral_taxa'],artificial_taxa='',cellular_taxa='',unknown_taxa=x['unknown_taxa'],near_best_taxids=''))
for name,data in [('comparison.tsv',rows),('eligible-nucleotide-alignments.tsv',details)]:
 with (d/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]) if data else ['query'],delimiter='\t');w.writeheader();w.writerows(data)
summary={'pairs':len(selection),'nr_complete':bool(a.nr),'groups':{},'transitions':{},'interpretation':'Sample-balanced secondary panel: one protein-only candidate beyond both NT24 and NT20 per sample, plus32 shared controls. Every mate is identical to its frozen NR input. Nucleotide classifications use top10 distinct score values and90% near-best per mate, then union mate taxa. Protein uses the unchanged NR summary. These are measured competing assignments, not truth labels or read-weighted accuracy.'}
tasks=sorted({x['task'] for x in rows});strata=sorted({x['stratum'] for x in rows});lookup={(x['task'],x['query_pair']):x for x in rows}
for task in tasks:
 for group in strata:
  selected=[x for x in rows if x['task']==task and x['stratum']==group]
  summary['groups'][task+'/'+group]=dict(pairs=len(selected),status=dict(collections.Counter(x['status'] for x in selected)),taxonomic_relation=dict(collections.Counter(x['relation_to_viral_reference'] for x in selected)))
for left,right in [('megablast','blastn')]+([('megablast','protein_NR'),('blastn','protein_NR')] if a.nr else []):
 for group in strata:
  selected=[q for q,x in selection.items() if x['stratum']==group]
  summary['transitions'][left+' -> '+right+'/'+group]=dict(pairs=len(selected),status=dict(collections.Counter(lookup[left,q]['status']+' -> '+lookup[right,q]['status'] for q in selected)),taxid_disagreements=sum(lookup[left,q]['taxid_lca']!=lookup[right,q]['taxid_lca'] for q in selected))
(d/'comparison.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

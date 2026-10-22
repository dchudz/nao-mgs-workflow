"""Retain every nucleotide rescue and describe strongest-alignment strength without reclassifying."""
import argparse,collections,csv,json,re,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('panel');a=p.parse_args();r=Path(a.panel)
by={(x['task'],x['query_pair']):x for x in csv.DictReader((r/'comparison.tsv').open(),delimiter='\t')}
selection={x['query_pair']:x for x in csv.DictReader((r/'selection.tsv').open(),delimiter='\t')}
ids={q for q in selection if selection[q]['stratum']=='added_beyond_NT24_and_NT20' and by['megablast',q]['status']=='no_alignment' and by['blastn',q]['status']=='viral_only'}
best={}
for x in csv.DictReader((r/'eligible-nucleotide-alignments.tsv').open(),delimiter='\t'):
 q=x['query'].rsplit('/',1)[0]
 if x['task']=='blastn' and q in ids and (q not in best or float(x['bitscore'])>float(best[q]['bitscore'])):best[q]=x
assert best.keys()==ids
members=collections.defaultdict(list)
for x in csv.DictReader((r/'original-members.tsv').open(),delimiter='\t'):members[x['query_pair']].append(x['seq_id'])
assert set(members)==set(selection)
rows=[dict(query_pair=q,sample=selection[q]['sample'],original_read_IDs=';'.join(members[q]),viral_reference_taxid=selection[q]['viral_reference_taxid'],blastn_taxid_lca=by['blastn',q]['taxid_lca'],relation_to_viral_reference=by['blastn',q]['relation_to_viral_reference'],aligned_mates=int(by['blastn',q]['aligned_mates']),**{k:best[q][k] for k in ['query','subject','bitscore','evalue','alignment_nt','identity_pct','coverage_pct']}) for q in sorted(ids)]
with (r/'rescued-protein-panel-pairs.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'queries':len(rows),'both_mates_aligned':sum(x['aligned_mates']==2 for x in rows),'relation_to_viral_reference':dict(collections.Counter(x['relation_to_viral_reference'] for x in rows)),'best_Evalue_above_1e_minus5':sum(float(x['evalue'])>1e-5 for x in rows),'interpretation':'Strongest reported eligible alignment per selected query pair. Exact-sequence query membership is distinct from independent searches. The1e-5 strength comparison uses the protein reporting cutoff descriptively and does not change nucleotide classifications. No biological truth or infection inferred.'}
for field,converter in [('alignment_nt',int),('identity_pct',float)]:
 values=[converter(x[field]) for x in rows]
 for stat,func in [('min',min),('median',statistics.median),('max',max)]:summary['best_'+field+'_'+stat]=func(values) if values else None
summary['largest_best_Evalue']=max((float(x['evalue']) for x in rows),default=None)
(r/'rescued-protein-panel-pairs.json').write_text(json.dumps(summary,indent=2))
resources=[]
for task in ['megablast','blastn']:
 x=json.loads((r/f'{task}.resources.json').read_text());text=(r/f'{task}.time').read_text()
 def number(label):return float(re.search(re.escape(label)+r':\s*([\d.]+)',text)[1])
 assert number('Exit status')==0
 resources.append(dict(task=task,queries=464,wall_seconds=x['wall_seconds'],cpu_seconds=number('User time (seconds)')+number('System time (seconds)'),peak_rss_kib=int(number('Maximum resident set size (kbytes)')),allocated_whole_worker_cost_usd=x['wall_seconds']/3600*(3.3345+(2048*.08+3000*.005+875*.04)/730)))
(r/'resource-summary.json').write_text(json.dumps({'runs':resources,'interpretation':'Shared primary worker and warm reference cache. Whole-worker allocation is not an incremental invoice charge. Timing includes all464panelmates, not just59rescues.'},indent=2));print(json.dumps(summary));print(json.dumps(resources))

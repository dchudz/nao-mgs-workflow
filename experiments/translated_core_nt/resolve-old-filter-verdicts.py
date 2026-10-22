import collections,csv,hashlib,json,shutil,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment');src=r/'cohort-storage/protein-nt-panel';out=r/'cohort-storage/translated-core-nt/old-panel-filter-audit/production-order';out.mkdir(exist_ok=False)
for name in ['selection.tsv','queries.fna']:shutil.copy2(src/name,out/name)
for method in ['megablast','blastn']:
 best={}
 for row in csv.reader((src/(method+'.tsv')).open(),delimiter='\t'):
  key=tuple(row[:2]);old=best.get(key)
  if old is None or (float(row[3]),int(row[5]))>(float(old[3]),int(old[5])):best[key]=row
 with (out/(method+'.tsv')).open('w') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(best.values())
cmd=['python3',str(r/'summarize-protein-nt-panel.py'),str(out),'--nodes',str(r/'reference/results/taxonomy-nodes.dmp'),'--selection-note','Filtering-order audit only: original unchanged searches, one best HSP per query/subject before the original rank filter.']
subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
before={(x['task'],x['query_pair']):x for x in csv.DictReader((src/'comparison.tsv').open(),delimiter='\t') if x['task']!='protein_NR'};after={(x['task'],x['query_pair']):x for x in csv.DictReader((out/'comparison.tsv').open(),delimiter='\t')};assert before.keys()==after.keys()
fields=['status','taxid_lca','relation_to_viral_reference'];changed=[]
for key,left in before.items():
 right=after[key]
 if any(left[k]!=right[k] for k in fields):changed.append({'method':key[0],'pair':key[1],'before':{k:left[k] for k in fields},'after':{k:right[k] for k in fields}})
result={'original_classification_sha256':hashlib.sha256((src/'comparison.tsv').read_bytes()).hexdigest(),'production_order_classification_sha256':hashlib.sha256((out/'comparison.tsv').read_bytes()).hexdigest(),'original_summary_source_sha256':hashlib.sha256((r/'summarize-protein-nt-panel.py').read_bytes()).hexdigest(),'changed_pair_verdicts':changed,'status_changes':sum(x['before']['status']!=x['after']['status'] for x in changed),'taxid_changes':sum(x['before']['taxid_lca']!=x['after']['taxid_lca'] for x in changed),'limits':['Only filtering order changes; same raw alignments, query bytes and taxonomy.','Historical reports and evidence are preserved; this is a labeled supplementary reanalysis.']};(out/'verdict-comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

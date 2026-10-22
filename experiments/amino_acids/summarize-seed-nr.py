"""Classify the same individual-mate sequences used in the nucleotide diagnostic."""
import collections,csv,json,sys
from pathlib import Path
sys.path.insert(0,'/home/ec2-user/experiment')
from importlib.machinery import SourceFileLoader
summarize=SourceFileLoader('nrsummary','/home/ec2-user/experiment/protein-nr-summary.py').load_module().summarize
r=Path('/home/ec2-user/experiment');d=Path(sys.argv[1]) if len(sys.argv)>1 else r/'protein/combined-validation';parents={}
for line in open(r/'protein/current-nr/nodes.dmp'):
 a=line.split('|');parents[int(a[0])]=int(a[1])
controls={'SEED_'+x['seq_id']:x for x in csv.DictReader(open(d/'seed-selection.tsv'),delimiter='\t')}
with open(d/'nr.tsv') as f,open(d/'seed-alignments.tsv','w') as selected:
 writer=csv.writer(selected,delimiter='\t')
 for row in csv.reader(f,delimiter='\t'):
  if row[0].startswith('SEED_'):writer.writerow(row)
with open(d/'seed-alignments.tsv') as f:
 rows=summarize(csv.reader(f,delimiter='\t'),parents,{q:x['original_taxid'] for q,x in controls.items()})
for x in rows:
 x['native_reference_taxid']=x.pop('viral_reference_taxid');x['sample']=controls[x['query_pair']]['sample'];x['seq_id']=controls[x['query_pair']]['seq_id']
with open(d/'seed-validation.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'diagnostic_pairs':len(rows),'status':dict(collections.Counter(x['nr_status'] for x in rows)),'relation':dict(collections.Counter(x['taxonomic_relation'] for x in rows)),'interpretation':'Same available mate sequences used in megablast and blastn diagnostic. NR protein coverage differs from core_nt nucleotide coverage; these are not known-truth labels.'}
(d/'seed-validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

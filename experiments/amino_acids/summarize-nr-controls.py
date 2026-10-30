"""Label baseline controls separately from the experimental protein candidates."""
import collections,csv,json,sys
from pathlib import Path
sys.path.insert(0,'/home/ec2-user/experiment')
from nr_evidence_inputs import alignment_input
from importlib.machinery import SourceFileLoader
summarize=SourceFileLoader('nrsummary','/home/ec2-user/experiment/protein-nr-summary.py').load_module().summarize
r=Path('/home/ec2-user/experiment');d=Path(sys.argv[1]) if len(sys.argv)>1 else r/'protein/combined-validation';parents={}
for line in open(r/'protein/current-nr/nodes.dmp'):
 a=line.split('|');parents[int(a[0])]=int(a[1])
controls={x['query_pair']:x for x in csv.DictReader(open(d/'native-controls.tsv'),delimiter='\t')}
with open(alignment_input(d,'native_controls')) as f:
 rows=summarize((x for x in csv.reader(f,delimiter='\t') if x[0].startswith('NR_CONTROL_')),parents,{q:x['native_taxid'] for q,x in controls.items()})
for x in rows:
 x['native_reference_taxid']=x.pop('viral_reference_taxid');x['sample']=controls[x['query_pair']]['sample'];x['example_seq_id']=controls[x['query_pair']]['example_seq_id']
with open(d/'native-control-validation.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['query_pair'],delimiter='\t');w.writeheader();w.writerows(rows)
summary={'control_pairs':len(rows),'status':dict(collections.Counter(x['nr_status'] for x in rows)),'relation':dict(collections.Counter(x['taxonomic_relation'] for x in rows)),'interpretation':'Secondary matched search of baseline nucleotide calls. Environmental calls are not known-truth positives. Selection is balanced by sample, not by viral species or read abundance.'}
(d/'native-control-validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

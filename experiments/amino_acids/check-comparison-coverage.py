"""Selected reads remain represented even when every method has no alignment."""
import csv,json,subprocess,tempfile
from pathlib import Path
r=Path(__file__).resolve().parent

def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

with tempfile.TemporaryDirectory() as tmp:
 d=Path(tmp);base=[dict(seq_id=q,sample='sample',AA12_added='True',NT20_added='False',baseline_control='False',original_taxid='10239') for q in ['q1','q2']]
 write(d/'selection.tsv',[dict(x,available_mates=2) for x in base])
 nt=[dict(x,task=task,status='no_alignment',lca=0,aligned_mates=0) for x in base for task in ['megablast','blastn']]
 nr=[dict(query_pair='SEED_'+x['seq_id'],nr_status='no_alignment',nr_taxid_lca=0,aligned_mates=0) for x in base]
 cmd=['python3',str(r/'compare-validation-methods.py'),'--nucleotide',str(d/'nt.tsv'),'--protein',str(d/'nr.tsv'),'--output',str(d/'out')]
 write(d/'nt.tsv',nt);write(d/'nr.tsv',nr)
 subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL);s=json.loads((d/'out/summary.json').read_text());assert s['pairs']==2 and s['available_query_sequences']==4
 # Removing the same selected read from all three methods must not pass the join.
 write(d/'nt.tsv',[x for x in nt if x['seq_id']=='q1']);write(d/'nr.tsv',[nr[0]])
 z=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);assert z.returncode and 'omits a selected original read pair' in z.stderr
 panel=d/'panel';panel.mkdir();selection=[dict(query_pair=q,sample='sample',stratum='control',viral_reference_taxid=10239) for q in ['q1','q2']];write(panel/'selection.tsv',selection)
 (panel/'queries.fna').write_text(''.join('>'+q+'/'+str(m)+'\nACGT\n' for q in ['q1','q2'] for m in [1,2]))
 for task in ['megablast','blastn']:(panel/(task+'.tsv')).write_text('')
 (d/'nodes.dmp').write_text('1\t|\t1\t|\tno rank\t|\n10239\t|\t1\t|\tsuperkingdom\t|\n')
 cmd=['python3',str(r/'summarize-protein-nt-panel.py'),str(panel),'--nodes',str(d/'nodes.dmp')]
 subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL);s=json.loads((panel/'comparison.json').read_text());assert s['pairs']==2 and all(z['status']=={'no_alignment':2} for z in s['groups'].values())
 write(panel/'selection.tsv',selection+[selection[0]]);z=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);assert z.returncode and 'Duplicate panel query' in z.stderr
 write(panel/'selection.tsv',selection);(panel/'megablast.tsv').write_text('q1/3\tsubject\t10239\t100\t1e-20\t30\t100\t100\n');z=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);assert z.returncode and 'assert q in set_query_ids' in z.stderr
print('Unaligned reads retained; jointly omitted reads, duplicate selections and unexpected mate IDs rejected.')

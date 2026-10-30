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
with tempfile.TemporaryDirectory() as tmp:
 d=Path(tmp)
 write(d/'evidence.tsv',[dict(query_pair='q1',exact_pair_also_native='True',original_pairs=2,added_original_pairs=0,nr_status='mixed_viral_nonviral'),dict(query_pair='q2',exact_pair_also_native='False',original_pairs=3,added_original_pairs=3,nr_status='viral_only')])
 write(d/'screens.tsv',[dict(query_pair='q1',original_members=2,nt_retained_original_members=2,nt_screen_status='all_retained',nt20_retained_original_members=2,nt20_screen_status='all_retained'),dict(query_pair='q2',original_members=3,nt_retained_original_members=0,nt_screen_status='none_retained',nt20_retained_original_members=1,nt20_screen_status='some_retained')])
 cmd=['python3',str(r/'analyze-protein-overlap.py'),str(d/'evidence.tsv'),str(d/'summary.json'),'--screen-membership',str(d/'screens.tsv')]
 result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);assert result.returncode and 'Specify the native arm' in result.stderr
 for arm,status,retained in [('nt','none_retained',0),('nt20','some_retained',1)]:
  subprocess.run(cmd+['--native-arm',arm],check=True,stdout=subprocess.DEVNULL);z=json.loads((d/'summary.json').read_text());group=z['groups']['viral_only'];category='new_exact_pair_all_original_members_absent_from_native'
  assert group['initial_screen/'+category+'/'+status]==1
  assert group['initial_screen_original_members/'+category+'/retained']==retained
  assert group['initial_screen_original_members/'+category+'/excluded']==3-retained
print('Coverage guards and arm-specific raw-screen boundary counts verified.')

"""Complete pair sequence identity, not read-ID multiplicity, determines exact-pair gain."""
import csv,gzip,json,subprocess,sys,tempfile
from pathlib import Path
src=Path(__file__).resolve().parent
union_source=src/'experiment-code/union_pairs.py'
if not union_source.exists():union_source=src/'union_pairs.py'
with tempfile.TemporaryDirectory(prefix='aa-comparison-') as tmp:
 r=Path(tmp);(r/'experiments/amino_acids').mkdir(parents=True);(r/'experiments/amino_acids/union_pairs.py').write_bytes(union_source.read_bytes())
 (r/'reference/results').mkdir(parents=True);(r/'reference/results/taxonomy-nodes.dmp').write_text('1 | 1 | no rank |\n10239 | 1 | superkingdom |\n100 | 10239 | species |\n')
 (r/'heldout-ranges.json').write_text(json.dumps([{'run_accession':'SAMPLE'}]));b=r/'cohort-storage/batches/b1';(b/'output/fastp').mkdir(parents=True);(b/'output/marked').mkdir();(b/'complete.json').write_text(json.dumps({'samples':['SAMPLE']}));p=r/'protein';p.mkdir()
 def table(path,rows):
  with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 table(b/'output/marked/x_SAMPLE_nt_duplicate_reads.tsv.gz',[{'seq_id':q,'aligner_taxid_lca':'100'} for q in ['old1','old2','lost']])
 with gzip.open(b/'output/fastp/SAMPLE_nt_fastp.fastq.gz','wt') as f:
  for q,s in [('old1','AAAA'),('old2','AAAA'),('lost','TTTT')]:
   for mate in [1,2]:f.write(f'@{q}/{mate}\n{s}\n+\nIIII\n')
 table(p/'all.tsv.gz',[{'query_pair':q,'example_seq_id':seqid,'occurrences':'1','query_seq':s,'query_seq_rev':s} for q,seqid,s in [('SAMPLE_p0','newID','AAAA'),('SAMPLE_p1','added','CCCC'),('SAMPLE_p2','old1','GGGG')]])
 table(p/'all.members.tsv.gz',[{'query_pair':'SAMPLE_p0','seq_id':'newID'},{'query_pair':'SAMPLE_p1','seq_id':'added'},{'query_pair':'SAMPLE_p2','seq_id':'old1'}])
 table(p/'viral-pairs.tsv',[{'query_pair':q,'protein_taxid_lca':'100','aligned_mates':'2','maximum_identity_pct':'100','maximum_alignment_aa':'25'} for q in ['SAMPLE_p0','SAMPLE_p1','SAMPLE_p2']])
 table(p/'nr.tsv',[{'query_pair':q,'nr_status':status,'aligned_mates':'2','taxonomic_relation':'exact','nr_taxid_lca':'100'} for q,status in [('SAMPLE_p0','viral_only'),('SAMPLE_p1','nonviral_only'),('SAMPLE_p2','viral_only')]])
 script=(src/'analyze-protein200.py').read_text().replace("Path('/home/ec2-user/experiment')",'Path('+repr(str(r))+')');(r/'analysis.py').write_text(script)
 subprocess.run([sys.executable,str(r/'analysis.py'),'--protein',str(p),'--nr',str(p/'nr.tsv'),'--output',str(r/'out')],check=True,stdout=subprocess.DEVNULL)
 result=json.loads((r/'out/summary.json').read_text())['tiers'];c=result['viral_reference_candidate'];n=result['nr_viral_only']
 assert (c['native_original_pairs'],c['native_exact_pairs'],c['added_exact_pairs'],c['lost_exact_pairs'])==(3,2,2,1),c
 assert (n['protein_exact_pairs'],n['added_exact_pairs'],n['shared_exact_pairs'],n['lost_exact_pairs'],n['added_original_pairs'])==(2,1,1,1,1),n
 subprocess.run([sys.executable,str(src/'analyze-protein-overlap.py'),str(r/'out/pair-evidence.tsv'),str(r/'overlap.json')],check=True,stdout=subprocess.DEVNULL)
 overlap=json.loads((r/'overlap.json').read_text())['groups']['all']
 assert overlap['exact_pair_also_native']==1
 assert overlap['new_exact_pair_all_original_members_absent_from_native']==1
 assert overlap['new_exact_pair_has_native_original_member']==1
 print('Verified: duplicate baseline IDs count once; a new ID with the same complete pair is not an exact-pair gain; nonviral competition removes only its own candidate.')

"""Isolate an observed score-statistics difference; reduced reference is only a mechanism fixture."""
import csv,hashlib,importlib.util,json,subprocess
from pathlib import Path
spec=importlib.util.spec_from_file_location('experiment',Path(__file__).with_name('experiment.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
d=m.ROOT/'cohort-storage/translated-core-nt/translation-fixture/word4-sum-statistics';d.mkdir(exist_ok=False)
query='SRR36995843.1123277743/2';seq=m.fasta(m.ROOT/'cohort-storage/translated-core-nt/pilot24/queries.fna')[query];(d/'query.fna').write_text('>'+query+'\n'+seq+'\n')
def docker(args):return ['docker','run','--rm','--network','none','--cpus','1','--memory','2g','-v',str(m.ROOT)+':'+str(m.ROOT),m.IMAGE]+args
extract=docker(['blastdbcmd','-db',str(m.DB),'-entry','MT766366.1','-outfmt','%f']);reference=subprocess.check_output(extract,text=True);(d/'reference.fna').write_text(reference)
subprocess.run(docker(['makeblastdb','-in',str(d/'reference.fna'),'-dbtype','nucl','-parse_seqids','-out',str(d/'db')]),check=True,stdout=subprocess.DEVNULL)
fields=m.FIELDS+['score'];variants={};commands={}
for word in ['default','word4']:
 for mode in ['default','true','false']:
  name=word+'-sum-'+mode;extra=(['-word_size','4','-threshold','18'] if word=='word4' else [])+(['-sum_stats',mode] if mode!='default' else [])
  cmd=docker(['tblastx']+extra+['-db',str(d/'db'),'-query',str(d/'query.fna'),'-query_gencode','1','-db_gencode','1','-searchsp','10000000000000','-qcov_hsp_perc','30','-evalue','10','-max_hsps','5','-num_alignments','250','-num_threads','1','-outfmt','6 '+' '.join(fields)]);run=subprocess.run(cmd,capture_output=True,text=True,timeout=60);assert run.returncode==0,run.stderr;(d/(name+'.tsv')).write_text(run.stdout);(d/(name+'.stderr')).write_text(run.stderr);commands[name]=cmd;variants[name]=[dict(zip(fields,line.split('\t'))) for line in run.stdout.splitlines()]
keyfields=['qseqid','sseqid','qstart','qend','sstart','send','qframe','sframe','qseq','sseq'];comparisons={}
for mode in ['default','true','false']:
 left={tuple(x[k] for k in keyfields):x for x in variants['default-sum-'+mode]};right={tuple(x[k] for k in keyfields):x for x in variants['word4-sum-'+mode]};shared=left.keys()&right.keys();comparisons[mode]={'shared_HSPs':len(shared),'different_bitscores':sum(left[k]['bitscore']!=right[k]['bitscore'] for k in shared),'different_Evalues':sum(left[k]['evalue']!=right[k]['evalue'] for k in shared),'different_raw_scores':sum(left[k]['score']!=right[k]['score'] for k in shared)}
result={'created_utc':m.now(),'excluded_from_real_data':True,'query':query,'reference_entry':'MT766366.1','reference_extraction_command':extract,'query_sha256':m.sha(d/'query.fna'),'reference_sha256':m.sha(d/'reference.fna'),'source_sha256':m.sha(__file__),'fixed_effective_search_space':10000000000000,'commands':commands,'alignments':variants,'comparisons':comparisons,'limits':['One observed query/subject pair, reduced-reference mechanism fixture only.','The arbitrary fixed search space makes within-fixture comparisons controlled; these Evalues are not validation against comprehensive core_nt.','No production verdict, experimental threshold or running search is changed.','Any mechanism demonstrated here must not be generalized to every observed Evalue difference without further evidence.']};m.save(d/'characterization.json',result);print(json.dumps({'comparisons':comparisons,'rows':{k:len(v) for k,v in variants.items()}}))

"""Inspect two already-returned reference proteins; never reclassify reads."""
import csv,datetime,hashlib,json,subprocess,time
from pathlib import Path
r=Path('/home/ec2-user/experiment');root=r/'protein/seed-nr-forced-double-index-trace/output-audit-v2';out=root/'reference-label-inspection';out.mkdir(exist_ok=False)
image='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
container=['docker','run','--rm','--network','none','--cpus','1','--memory','4g','-v',str(r)+':'+str(r)+':ro',image]
def run(cmd,name,seconds):
    start=time.monotonic();result=subprocess.run(cmd,capture_output=True,text=True,timeout=seconds)
    (out/(name+'.out')).write_text(result.stdout);(out/(name+'.log')).write_text(result.stderr)
    receipt=dict(command=cmd,returncode=result.returncode,wall_seconds=time.monotonic()-start)
    (out/(name+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');assert result.returncode==0
    return result.stdout,receipt
viral,fetch=run(container+['blastdbcmd','-db',str(r/'protein/current-nr/nr'),'-entry','BDA36041.1','-outfmt','%f'],'viral-retrieval',120)
viral_file=out/'viral-BDA36041.1.faa';viral_file.write_text(viral)
human_file=root/'human-target-native-retrieval.faa';human=human_file.read_text()
for text in (viral,human):assert text.startswith('>') and text.count('>')==1
fields='qseqid sseqid pident length qlen slen qstart qend sstart send bitscore evalue qseq sseq'.split()
cmd=container+['blastp','-query',str(human_file),'-subject',str(viral_file),'-num_threads','1','-max_hsps','1','-outfmt','6 '+' '.join(fields)]
text,alignment=run(cmd,'reference-alignment',120)
rows=[dict(zip(fields,x)) for x in csv.reader(text.splitlines(),delimiter='\t')];assert all(len(x.split('\t'))==len(fields) for x in text.splitlines())
for x in rows:
    x['query_span_fraction']=(abs(int(x['qend'])-int(x['qstart']))+1)/int(x['qlen'])
    x['subject_span_fraction']=(abs(int(x['send'])-int(x['sstart']))+1)/int(x['slen'])
receipt=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),human_header=human.splitlines()[0],viral_header=viral.splitlines()[0],human_fasta_sha256=hashlib.sha256(human_file.read_bytes()).hexdigest(),viral_fasta_sha256=hashlib.sha256(viral_file.read_bytes()).hexdigest(),retrieval=fetch,alignment=alignment,rows=rows,limits=['Descriptive pairwise inspection of two already-returned NR targets; not comprehensive validation and not a new read verdict.','This small subject set changes E-value calibration. Its E-values must not be compared with comprehensive NR or core_nt searches.','Identity and alignment spans describe these two stored sequences, not the origin of the reads or which reference annotation is correct.','Original human and viral taxonomic labels, near-best competition, and all reported candidate counts remain unchanged.'])
(out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k not in ['rows','retrieval','alignment']},indent=2));print(json.dumps([{k:v for k,v in x.items() if k not in ['qseq','sseq']} for x in rows]))

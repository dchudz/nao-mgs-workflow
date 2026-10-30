"""Locate previously validated case reads on one existing viral reference per sample."""
import argparse,collections,csv,datetime,hashlib,json,subprocess,time
from pathlib import Path
from Bio.Seq import Seq
IMAGE='public.ecr.aws/q0n1c7g8/nao-mgs-workflow/blast@sha256:cb02427c58f3842e371a8dfea28d2a92e311406dbebec2c3a457e375079f2602'
FIELDS='qseqid sseqid bitscore evalue length pident qlen qstart qend sstart send qframe sframe qseq sseq'.split()
def fasta(p):
 d={};name=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];d[name]=''
  else:d[name]+=line.strip()
 return d
def extract(seq,start,end,frame):
 nt=seq[min(start,end)-1:max(start,end)]
 return str(Seq(nt).reverse_complement()) if frame<0 else nt
def regions(intervals):
 out=[]
 for lo,hi in sorted(intervals):
  if out and lo<=out[-1][1]:out[-1][1]=max(hi,out[-1][1])
  else:out.append([lo,hi])
 return out
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--commit',required=True);a=p.parse_args();protocol=json.loads((a.directory/'protocol.json').read_text());commands=[];summaries=[];allrows=[];start=time.monotonic()
for case in protocol['cases']:
 d=a.directory/case['sample'];q=d/'queries.fna';ref=d/'reference.fna';assert hashlib.sha256(q.read_bytes()).hexdigest()==case['queries_sha256'];assert hashlib.sha256(ref.read_bytes()).hexdigest()==case['reference_sha256']
 out=d/'locus-alignments.tsv';assert not out.exists()
 prefix=['docker','run','--rm','--network','none','--cpus','1','--memory','4g','-v',str(a.directory)+':'+str(a.directory),IMAGE]
 make=prefix+['makeblastdb','-in',str(ref),'-dbtype','nucl','-parse_seqids','-out',str(d/'reference-db')]
 command=prefix+['tblastx','-query',str(q),'-db',str(d/'reference-db'),'-query_gencode','1','-db_gencode','1','-qcov_hsp_perc','30','-evalue','10','-num_threads','1','-max_hsps','5','-outfmt','6 '+' '.join(FIELDS),'-out',str(out)]
 with (d/'search.log').open('w') as log:
  subprocess.run(make,stdout=log,stderr=log,check=True,timeout=120);t=time.monotonic();subprocess.run(command,stdout=log,stderr=log,check=True,timeout=120);wall=time.monotonic()-t
 commands.extend([make,command]);queries=fasta(q);reference=next(iter(fasta(ref).values()));best={}
 for row in csv.DictReader(out.open(),fieldnames=FIELDS,delimiter='\t'):
  if int(row['length'])<25:continue
  key=row['qseqid'];rank=lambda x:(-float(x['bitscore']),-int(x['length']),int(x['sstart']),int(x['send']),int(x['qstart']),int(x['qend']))
  if key not in best or rank(row)<rank(best[key]):best[key]=row
 rows=[]
 for query,row in sorted(best.items()):
  qnt=extract(queries[query],int(row['qstart']),int(row['qend']),int(row['qframe']));snt=extract(reference,int(row['sstart']),int(row['send']),int(row['sframe']));qaa=str(Seq(qnt).translate());saa=str(Seq(snt).translate());assert qaa.upper()==row['qseq'].upper() and saa.upper()==row['sseq'].upper();assert len(qnt)==len(snt)==int(row['length'])*3
  row.update(sample=case['sample'],reference=case['reference'],start_0based=min(int(row['sstart']),int(row['send']))-1,end_exclusive=max(int(row['sstart']),int(row['send'])),nt_identity_pct=100*sum(x==y for x,y in zip(qnt,snt))/len(qnt),synonymous_changed_codons=sum(x==y and x in 'ACDEFGHIKLMNPQRSTVWY' and qnt[i*3:i*3+3]!=snt[i*3:i*3+3] for i,(x,y) in enumerate(zip(qaa,saa))));rows.append(row)
 merged=regions([(x['start_0based'],x['end_exclusive']) for x in rows]);summary=dict(case,search_wall_seconds=wall,aligned_mates_at_least25aa=len(rows),audited_pairs=len({q.rsplit('/',1)[0] for q in queries}),aligned_pairs=len({x['qseqid'].rsplit('/',1)[0] for x in rows}),reference_intervals_0based_half_open=merged,reference_covered_nt=sum(hi-lo for lo,hi in merged),mates_with_synonymous_changes=sum(x['synonymous_changed_codons']>0 for x in rows));summaries.append(summary);allrows+=rows
with (a.directory/'best-locus-alignments.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':a.commit,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.monotonic()-start,'cases':summaries,'commands':commands,'limits':protocol['limits']+['Best eligible ungapped TBLASTX HSP per mate, at least25aminoacids and30percent query coverage. Every query retained in the denominator.','Intervals are overlaps on a single selected reference; neither a production duplicate check nor proof of independent molecules.','This audit leaves every comprehensive NR and core_nt validation result unchanged.']};(a.directory/'locus-summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(summaries))

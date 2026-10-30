"""Record the actual full-cohort validation queries and search traces, including command provenance."""
import csv,gzip,hashlib,json,re
from pathlib import Path
r=Path('/home/ec2-user/experiment');d=r/'cohort-storage/validation-cohort'
assert json.loads((d/'complete.json').read_text())['completed_samples']==1000
units={'ms':.001,'s':1,'m':60,'h':3600,'d':86400}
def seconds(s):
 if s=='-':return None
 parts=re.findall(r'([\d.]+)(ms|s|m|h|d)',s);assert parts,s
 return sum(float(n)*units[u] for n,u in parts)
rows=[]
for x in csv.DictReader((d/'output/trace.tsv').open(),delimiter='\t'):
 m=re.fullmatch(r'VALIDATE_VIRAL_ASSIGNMENTS:BLAST_FASTA:BLAST \(id=(nt|nt22|nt20|aa12|union)\)',x['name'])
 if not m:continue
 arm=m[1];assert x['status']=='COMPLETED' and x['exit']=='0'
 matches=list((d/'work').glob(x['hash']+'*'));assert len(matches)==1
 work=matches[0];query=work/(arm+'_sampled_reads.fasta.gz');assert query.exists()
 h=hashlib.sha256();query_count=0
 with gzip.open(query,'rb') as f:
  for line in f:h.update(line);query_count+=line.startswith(b'>')
 elapsed=seconds(x['realtime']);cpu=None if x['%cpu']=='-' or elapsed is None else elapsed*float(x['%cpu'].rstrip('%'))/100
 rows.append({'arm':arm,'queries':query_count,'uncompressed_query_sha256':h.hexdigest(),'task_realtime_seconds':elapsed,'task_CPU_seconds_estimate':cpu,'task_peak_rss':x['peak_rss'],'command_script':(work/'.command.sh').read_text(),'work_hash':x['hash'],'allocated_worker_cost_usd':None if elapsed is None else elapsed/3600*(3.3345+(2048*.08+3000*.005+875*.04)/730)})
assert len(rows)==5 and {x['arm'] for x in rows}=={'nt','nt22','nt20','aa12','union'}
order=['nt','nt22','nt20','aa12','union'];rows.sort(key=lambda x:order.index(x['arm']))
summary={'completed_samples':1000,'runs':rows,'interpretation':['Original joined/merged-mate production validation with20selectedexemplars per species per pooled arm.','Realtimes and CPU estimates come from completed Nextflow search tasks; they include task startup and the cached database availability check.','The default blastn task is megablast. Allfivearms use16threads and the same comprehensive core_nt.','Whole-worker allocations include provisioned EBS but are not incremental charges and must not be summed with overlapping-job allocations.','The failed orchestration attempt is preserved separately; completed search tasks were first executed after the explicit-session resume.']}
(d/'resource-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps([{k:v for k,v in x.items() if k!='command_script'} for x in rows]))

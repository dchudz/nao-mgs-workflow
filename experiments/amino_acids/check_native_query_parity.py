from pathlib import Path
import csv,gzip,hashlib,json
r=Path('/home/ec2-user/experiment/e2e/SRR37076045');data={}
for arm in ['production','validation_control']:
 rows=list(csv.DictReader((r/(arm+'-downstream')/'trace.tsv').open(),delimiter='\t'))
 row=next(x for x in rows if x['process'].endswith(':BLAST'))
 d=Path(row['workdir']);p=d/f'{r.name}_in.fasta.gz'
 with gzip.open(p,'rb') as f:b=f.read()
 data[arm]={'sha256':hashlib.sha256(b).hexdigest(),'queries':sum(x.startswith(b'>') for x in b.splitlines()),'input_name':p.name,'blast_trace':{k:row[k] for k in ['realtime','%cpu','peak_rss','read_bytes','write_bytes']}}
assert len(data)==2 and data['production']['sha256']==data['validation_control']['sha256']
record={'verified':True,'sample':r.name,'arms':data,'timing_limit':'Task times include concurrent AA screening and differing database cache state; CPU and I/O measurements are retained alongside wall time.'}
(r/'validation-query-parity.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))

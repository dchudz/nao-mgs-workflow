"""Record cumulative live CPU counters so an incomplete trial retains a measured lower bound."""
import datetime,json,os,subprocess
from pathlib import Path
r=Path('/home/ec2-user/experiment/cohort-storage/blast-threading');rows=[]
for i in [1,2,3]:
 name=f'aa-blast-threading-{i}';p=subprocess.run(['docker','top',name,'-eo','pid,comm'],capture_output=True,text=True)
 if p.returncode:continue
 for line in p.stdout.splitlines()[1:]:
  pid,comm=line.split()
  if comm!='blastn':continue
  root=Path('/proc')/pid
  try:
   fields=(root/'stat').read_text().rsplit(')',1)[1].split();ticks=os.sysconf('SC_CLK_TCK');status={line.split(':',1)[0]:line.split(':',1)[1].strip() for line in (root/'status').read_text().splitlines()}
  except FileNotFoundError:continue
  rows.append({'container':name,'PID':int(pid),'CPU_seconds_observed':(int(fields[11])+int(fields[12]))/ticks,'user_CPU_seconds':int(fields[11])/ticks,'system_CPU_seconds':int(fields[12])/ticks,'VmHWM':status.get('VmHWM'),'VmRSS':status.get('VmRSS'),'scope':'Cumulative process CPU counter while running, including process threads; lower bound if no completed GNU time measurement is available.'})
snapshot={'observed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'processes':rows}
with (r/'live-resource-observations.jsonl').open('a') as f:f.write(json.dumps(snapshot)+'\n')
print(json.dumps(snapshot))

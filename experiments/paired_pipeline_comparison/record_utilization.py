"""Save existing CloudWatch measurements; do not alter worker configuration."""
import argparse
import datetime as dt
import json
import statistics
import subprocess
from pathlib import Path

INSTANCES={'spot2resume':'i-0638e9b2262f8f52c','primary':'i-02d118993a6f4989a','spot1':'i-0f5065e3b4f67cfbc','spot2':'i-0ee9da3f09d197be6'}
VOLUMES={'resumeboot':'vol-05bffc9a8717b6229','primaryroot':'vol-0272ade891333af83','spot1root':'vol-004c4ffab08c46dbd','spot2root':'vol-05ca1445d265f5325','spot1scratch':'vol-06f002b044d0aeb35','spot2scratch':'vol-03f6a458f0812c7b4'}
AWS=['aws','--profile','pipeline-codex','--region','us-east-1']

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);a=parser.parse_args()
    end=dt.datetime.now(dt.timezone.utc);start='2026-10-24T05:00:00Z';queries=[]
    for worker,instance in INSTANCES.items():
        queries.append({'Id':worker+'cpu','Label':worker+' CPU percent','MetricStat':{'Metric':{'Namespace':'AWS/EC2','MetricName':'CPUUtilization','Dimensions':[{'Name':'InstanceId','Value':instance}]},'Period':300,'Stat':'Average'},'ReturnData':True})
    for name,volume in VOLUMES.items():
        for metric,stat,suffix in [('VolumeReadBytes','Sum','read'),('VolumeWriteBytes','Sum','write'),('VolumeQueueLength','Average','queue')]:
            queries.append({'Id':name+suffix,'Label':name+' '+metric,'MetricStat':{'Metric':{'Namespace':'AWS/EBS','MetricName':metric,'Dimensions':[{'Name':'VolumeId','Value':volume}]},'Period':300,'Stat':stat},'ReturnData':True})
    a.out.parent.mkdir(parents=True,exist_ok=True);queryfile=a.out.with_suffix('.queries.json');queryfile.write_text(json.dumps(queries,indent=2)+'\n')
    command=AWS+['cloudwatch','get-metric-data','--metric-data-queries','file://'+str(queryfile.resolve()),'--start-time',start,'--end-time',end.isoformat(),'--scan-by','TimestampAscending','--max-datapoints','100800']
    data=json.loads(subprocess.check_output(command,text=True));assert not data.get('NextToken')
    summary=[]
    for row in data['MetricDataResults']:
        assert row['StatusCode']=='Complete';values=row['Values'];assert len(values)==len(row['Timestamps'])
        if values:summary.append({'metric_id':row['Id'],'label':row['Label'],'points':len(values),'mean_of_reported_periods':statistics.mean(values),'maximum_reported_period':max(values)})
    record={'start_utc':start,'end_utc':end.isoformat(),'period_seconds':300,'queries':queries,'response':data,'summary':summary,'limits':['Missing periods are not zero. Means cover reported periods only, including preparation and idle periods.','Five-minute intervals can span several short pipeline stages; these are not per-stage CPU measurements.','CPU percentage is the instance-wide fraction across 48 vCPUs, not the CPU allocation requested by a task.','Storage metrics record device traffic and queueing, not a causal attribution of stage runtime. No memory metric is available here.']}
    a.out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'metrics':len(summary),'points':sum(r['points'] for r in summary),'cpu':[r for r in summary if r['metric_id'].endswith('cpu')]}))

if __name__=='__main__':main()

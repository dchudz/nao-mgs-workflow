"""Measure AA child completion separately from a paused parent orchestrator."""
import argparse,datetime,json,os,re,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('sample');p.add_argument('pid',type=int);a=p.parse_args()
R=Path('/home/ec2-user/experiment');D=R/'e2e'/a.sample;proc=Path('/proc')/str(a.pid)
cmd=proc.joinpath('cmdline').read_bytes().split(b'\0');assert cmd[0].endswith(b'/screen_amino_acids')
assert any(a.sample.encode() in x for x in cmd)
stat=proc.joinpath('stat').read_text().rsplit(') ',1)[1].split();hz=os.sysconf('SC_CLK_TCK')
start_ticks=int(stat[19]);boot=time.clock_gettime(time.CLOCK_BOOTTIME)
started=datetime.datetime.now(datetime.timezone.utc)-datetime.timedelta(seconds=boot-start_ticks/hz)
record={'pid':a.pid,'command':[x.decode() for x in cmd if x],'started_utc_from_kernel_clock':started.isoformat(),
 'clock_resolution_seconds':1/hz,'scientific_source_commit':(R/'source-commit.txt').read_text().strip(),
 'reason':'Primary orchestrator may be paused while independent native controls run. Its stage receipt can include waiting after this child has finished.'}
expected=json.loads((R/'data'/a.sample/'verified.json').read_text())['spot_count']
log=D/'aa-screen.log';last_cpu=None;final_cpu=False
while True:
    try:
        fields=proc.joinpath('stat').read_text().rsplit(') ',1)[1].split()
        last_cpu=(int(fields[11])+int(fields[12]))/hz;final_cpu=fields[0]=='Z'
    except FileNotFoundError:fields=None
    text=log.read_text();match=re.search(r'matched_pairs=(\d+) total_pairs=(\d+) amino_k=(\d+)',text)
    if match:
        matched,total,k=map(int,match.groups());assert total==expected and k==12
        finished=datetime.datetime.fromtimestamp(log.stat().st_mtime,datetime.timezone.utc)
        record.update(verified=True,matched_pairs=matched,total_pairs=total,amino_k=k,
            final_counts_written_utc=finished.isoformat(),wall_seconds_until_final_counts=(finished-started).total_seconds(),
            observed_CPU_seconds=last_cpu,CPU_is_final=final_cpu,
            CPU_note='Final process CPU if a zombie was observed; otherwise a lower bound from the latest live observation.')
        break
    if fields is None or fields[0]=='Z':raise RuntimeError('AA screen exited without successful final counts')
    time.sleep(5)
partial=D/'aa-screen-execution.partial';partial.write_text(json.dumps(record,indent=2)+'\n');partial.replace(D/'aa-screen-execution.json')
print(json.dumps(record),flush=True)

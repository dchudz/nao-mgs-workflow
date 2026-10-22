"""Foreground read-only wait; return when a new result needs scientific review."""
import concurrent.futures, datetime, json, shlex, subprocess, time
from pathlib import Path
r=Path(__file__).resolve().parent
remote='''import json,shutil,sys,re
from pathlib import Path
r=Path('/home/ec2-user/experiment');h=sys.argv[1];o={'host':h,'free_GiB':round(shutil.disk_usage(r).free/2**30,1)}
if h=='aa-screen':
 d=r/'cohort/heldout';o.update(screened=len(list(d.glob('*/complete.json'))),failures=len(list(d.glob('*/failure.json'))))
 o['protein_repeat_gate_done']=(r/'protein-repeat-gate/complete.json').exists()
 logs=['heldout-screen-rebalanced.log','protein-repeat-gate-rebalanced.log']
elif h=='aa-experiment':
 o['primary_screens']=len(list((r/'cohort/heldout').glob('*/complete.json')))
 o['primary_screen_failures']=len(list((r/'cohort/heldout').glob('*/failure.json')))
 o['native']=sum(len(json.loads(p.read_text())['samples']) for p in (r/'cohort-storage/batches').glob('*/complete.json'))
 o['repeat']=sum(len(json.loads(p.read_text())['samples']) for p in (r/'cohort-storage/repeat-confirmation/batches').glob('*/complete.json'))
 o['scaling_done']=(r/'cohort-storage/validation-scaling/complete.json').exists()
 o['protein_nt_panel_done']=(r/'cohort-storage/protein-nt-panel/complete.json').exists()
 o['scaling_runs']=len(list((r/'cohort-storage/validation-scaling').glob('run*/resources.json')))
 o['seed_done']=(r/'cohort-storage/seed-diagnostic200/blastn.resources.json').exists()
 o['validation_done']=(r/'cohort-storage/validation-cohort/complete.json').exists()
 o['producer_commit_done']=(r/'cohort-storage/producer-finished.json').exists()
 o['cohort_free_GiB']=round(shutil.disk_usage(r/'cohort-storage').free/2**30,1)
 logs=['heldout-primary-screen.log','validation-scaling.log','validation-scaling-resumed.log','heldout-consumer-recovery.log','heldout-seed-diagnostic.log','repeat-confirmation.log','cohort-validation-scheduled.log','await-screen-commit.log','protein-nt-panel.log']
else:
 o['protein_done']=(r/'cohort-storage/protein200/viral-pairs.tsv.json').exists()
 d=r/'protein/combined-validation';o['nr_started']=(d/'search-protocol.json').exists();o['nr_done']=(d/'complete.json').exists();o['nr_incomplete']=(d/'incomplete.json').exists()
 o['dedup_done']=(r/'protein/nr-dedup-benchmark/summary.json').exists()
 log=d/'nr-search.log'
 if log.exists():
  with log.open('rb') as f:f.seek(max(0,log.stat().st_size-8192));tail=f.read().decode(errors='replace')
  blocks=re.findall(r'reference block ([0-9]+)/([0-9]+)',tail)
  if blocks:o['nr_reference_block'],o['nr_reference_blocks']=map(int,blocks[-1])
 o['nr_reference_done']=(r/'protein/current-nr/download-complete.json').exists() if h=='aa-nr-fast' else False
 o['thread_parity_done']=(r/'protein/nr-thread-parity-on-demand/parity.json').exists()
 logs=['nr-on-demand.log','heldout-protein-preparation.log','nr-dedup-benchmark.log','cohort-storage/protein200/viral-search.log']
errors=[]
for name in logs:
 p=r/name
 if p.exists():
  with p.open('rb') as f:f.seek(max(0,p.stat().st_size-16384));tail=f.read().decode(errors='replace')
  if 'Traceback (most recent call last):' in tail or '\\nError:' in tail:errors.append(name)
o['error_logs']=errors
print(json.dumps(o))
'''
def read(host):
 try:
  return json.loads(subprocess.check_output(['ssh','-F',str(r/'ssh_config'),host,'python3 -c '+shlex.quote(remote)+' '+host],text=True,timeout=25))
 except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as error:
  return {'host':host,'free_GiB':1000,'connection_error':type(error).__name__}
def events(rows):
 result=set()
 for d in rows:
  for k,v in d.items():
   if k.endswith('_done') or k in ['nr_started','nr_incomplete']:
    if v:result.add(d['host']+'/'+k)
  if d.get('repeat',0)>0:result.add('first_repeat_batch')
  for count in range(1,d.get('scaling_runs',0)+1):result.add('scaling_run_'+str(count))
  if d.get('native',0)>=500:result.add('native_half_cohort')
  if d.get('screened',0)==800:result.add('screen_first800_done')
  if d.get('native',0)==1000:result.add('native_1000')
  if d.get('screened',0)==1000:result.add('screen_1000')
  if d.get('failures',0) or d.get('primary_screen_failures',0):result.add('screen_failure')
  if d.get('primary_screens',0)==200:result.add('primary_screen200_done')
  if d.get('connection_error'):result.add(d['host']+'_connection_error')
  for name in d.get('error_logs',[]):result.add(d['host']+'/job_error/'+name)
  if d.get('cohort_free_GiB',1000)<60:result.add('cohort_disk_attention')
  if d['free_GiB'] < (120 if d['host']=='aa-screen' else 60):result.add(d['host']+'_disk_attention')
 return result
state=r/'foreground-watch-seen.json';seen=set(json.loads(state.read_text())) if state.exists() else set()
start=time.monotonic();last_print=0
while True:
 with concurrent.futures.ThreadPoolExecutor(3) as ex:rows=list(ex.map(read,['aa-experiment','aa-screen','aa-nr']))
 found=events(rows)
 # Clear recovered alerts so a later recurrence wakes the foreground wait again.
 seen={key for key in seen if not any(tag in key for tag in ['job_error','connection_error','disk_attention','screen_failure']) or key in found}
 state.write_text(json.dumps(sorted(seen)))
 new=found-seen
 snapshot={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workers':rows,'new_events':sorted(new)}
 (r/'watch-status-latest.json').write_text(json.dumps(snapshot,indent=2))
 with (r/'progress-history.jsonl').open('a') as f:f.write(json.dumps(snapshot)+'\n')
 if new or time.monotonic()-last_print>300:
  print(json.dumps(snapshot),flush=True);last_print=time.monotonic()
 if new:
  state.write_text(json.dumps(sorted(seen|new)));break
 if time.monotonic()-start>900:break
 time.sleep(30)

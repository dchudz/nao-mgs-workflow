"""Foreground read-only wait; return when a new result needs scientific review."""
import concurrent.futures, datetime, json, shlex, subprocess, time
from pathlib import Path
r=Path(__file__).resolve().parent
remote='''import json,shutil,sys,re
from pathlib import Path
r=Path('/home/ec2-user/experiment');h=sys.argv[1];o={'host':h,'free_GiB':round(shutil.disk_usage(r).free/2**30,1)}
if h=='aa-screen':
 d=r/'cohort/heldout';o.update(screened=len(list(d.glob('*/complete.json'))),failures=len(list(d.glob('*/failure.json'))))
 o['screen_closure_done']=(r/'screen-closure/closure-backup-verification.json').exists()
 logs=['heldout-screen-rebalanced.log','screen-closure.log']
elif h=='aa-experiment':
 final_log=r/'finalize-cohort.log';o['native_analysis_done']=final_log.exists() and 'FULL1000_NATIVE_ANALYSES_COMPLETE' in final_log.read_text()
 o['protein_primary_analysis_done']=(r/'cohort-storage/protein-finalization/complete.json').exists()
 o['validation_started']=(r/'cohort-storage/validation-cohort/scope.json').exists()
 o['final_cohort_analyses_done']=(r/'cohort-storage/analysis/final-cohort-analyses.json').exists()
 o['protein_screen_boundaries_done']=(r/'protein-screen-boundaries/complete.json').exists()
 o['protein_repeat_gate_done']=(r/'protein-repeat-gate/complete.json').exists()
 o['primary_screens']=len(list((r/'cohort/heldout').glob('*/complete.json')))
 o['primary_screen_failures']=len(list((r/'cohort/heldout').glob('*/failure.json')))
 o['native']=sum(len(json.loads(p.read_text())['samples']) for p in (r/'cohort-storage/batches').glob('*/complete.json'))
 o['repeat']=sum(len(json.loads(p.read_text())['samples']) for p in (r/'cohort-storage/repeat-confirmation/batches').glob('*/complete.json'))
 o['scaling_done']=(r/'cohort-storage/validation-scaling/complete.json').exists()
 o['protein_nt_panel_done']=(r/'cohort-storage/protein-nt-panel/complete.json').exists()
 o['protein_nt_panel_runs']=len(list((r/'cohort-storage/protein-nt-panel').glob('*.resources.json')))
 o['header_diagnostic_done']=(r/'cohort-storage/preparation-discrepancies/nt-validation/complete.json').exists()
 o['blast_threading_done']=(r/'cohort-storage/blast-threading/complete.json').exists()
 o['blast_threading_runs']=len(list((r/'cohort-storage/blast-threading').glob('run*/resources.json')))
 o['scaling_runs']=len(list((r/'cohort-storage/validation-scaling').glob('run*/resources.json')))
 o['seed_done']=(r/'cohort-storage/seed-diagnostic200/blastn.resources.json').exists()
 o['validation_done']=(r/'cohort-storage/validation-cohort/complete.json').exists()
 o['producer_commit_done']=(r/'cohort-storage/producer-finished.json').exists()
 o['cohort_free_GiB']=round(shutil.disk_usage(r/'cohort-storage').free/2**30,1)
 logs=['heldout-primary-screen.log','validation-scaling.log','validation-scaling-resumed.log','heldout-consumer-recovery.log','heldout-seed-diagnostic.log','repeat-confirmation.log',('cohort-validation-resumed.log' if (r/'cohort-validation-resumed.log').exists() else 'cohort-validation-scheduled.log'),'await-screen-commit.log','protein-nt-panel.log','blast-threading.log','header-discrepancy-validation.log','protein-repeat-gate-primary.log','protein-screen-boundaries.log','finalize-cohort.log','finalize-protein-primary.log']
else:
 o['protein_done']=(r/'cohort-storage/protein200/viral-pairs.tsv.json').exists()
 d=r/'protein/combined-validation';protocol=d/'search-protocol.json';o['nr_started']=protocol.exists() and json.loads(protocol.read_text()).get('recovery_from')=='nr24-oom-20261022T034211';o['nr_done']=(d/'complete.json').exists();o['nr_incomplete']=(d/'incomplete.json').exists()
 o['dedup_done']=(r/'protein/nr-dedup-benchmark/summary.json').exists()
 log=d/'nr-search.log'
 if log.exists() and o['nr_started']:
  with log.open('rb') as f:f.seek(max(0,log.stat().st_size-8192));tail=f.read().decode(errors='replace')
  blocks=re.findall(r'reference block ([0-9]+)/([0-9]+)',tail)
  if blocks:o['nr_reference_block'],o['nr_reference_blocks']=map(int,blocks[-1])
 o['nr_reference_done']=(r/'protein/current-nr/download-complete.json').exists() if h=='aa-nr-fast' else False
 metrics=d/'memory-metrics.jsonl'
 if metrics.exists():
  with metrics.open('rb') as f:f.seek(max(0,metrics.stat().st_size-1024));last=f.read().decode().strip().splitlines()[-1]
  m=json.loads(last);o['nr_memory_available_GiB']=round(m['MemAvailable']/1024**2,2);o['nr_swap_used_GiB']=round((m['SwapTotal']-m['SwapFree'])/1024**2,2);o['nr_process_RSS_GiB']=round(m['diamond_RSS_KiB']/1024**2,2) if m['diamond_RSS_KiB'] is not None else None
 o['nr_primary_analysis_done']=(r/'nr-primary-analysis-complete.json').exists()
 o['nr_raw_backup_done']=(d/'search-raw-complete.json').exists()
 o['standalone_seed_nr_done']=(r/'protein/seed-nr-standalone/complete.json').exists()
 o['standalone_seed_nr_started']=(r/'protein/seed-nr-standalone/search-protocol.json').exists()
 o['standalone_seed_nr_incomplete']=(r/'protein/seed-nr-standalone/incomplete.json').exists()
 o['thread_parity_done']=(r/'protein/nr-thread-parity-on-demand/parity.json').exists()
 o['nr_memory_restart_parity_done']=(r/'protein/nr-thread-parity-memory-restart/parity.json').exists()
 logs=['finalize-protein-nr.log','nr-restart256.log','standalone-seed-nr.log','heldout-protein-preparation.log','nr-dedup-benchmark.log','cohort-storage/protein200/viral-search.log']
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
 expected=r/'expected-worker-states.json'
 state=json.loads(expected.read_text()).get(host) if expected.exists() else None
 if state:return {'host':host,'expected_state':state,'free_GiB':None}
 try:
  return json.loads(subprocess.check_output(['ssh','-F',str(r/'ssh_config'),host,'python3 -c '+shlex.quote(remote)+' '+host],text=True,timeout=25))
 except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as error:
  return {'host':host,'free_GiB':1000,'connection_error':type(error).__name__}
def events(rows):
 result=set()
 for d in rows:
  if d.get('expected_state'):continue
  for k,v in d.items():
   if k.endswith(('_done','_started','_incomplete')):
    if v:result.add(d['host']+'/'+k)
  if d.get('repeat',0)>0:result.add('first_repeat_batch')
  for count in range(1,d.get('scaling_runs',0)+1):result.add('scaling_run_'+str(count))
  for count in range(1,d.get('protein_nt_panel_runs',0)+1):result.add('protein_nt_panel_run_'+str(count))
  for count in range(1,d.get('blast_threading_runs',0)+1):result.add('blast_threading_run_'+str(count))
  if d.get('native',0)>=500:result.add('native_half_cohort')
  if d.get('screened',0)==800:result.add('screen_first800_done')
  if d.get('native',0)==1000:result.add('native_1000')
  if d.get('screened',0)==1000:result.add('screen_1000')
  if d.get('failures',0) or d.get('primary_screen_failures',0):result.add('screen_failure')
  if d.get('primary_screens',0)==200:result.add('primary_screen200_done')
  if d.get('connection_error'):result.add(d['host']+'_connection_error')
  for name in d.get('error_logs',[]):result.add(d['host']+'/job_error/'+name)
  if d.get('nr_memory_available_GiB',1000)<12 or d.get('nr_swap_used_GiB',0)>16:result.add('NR_memory_attention')
  if d.get('cohort_free_GiB',1000)<60:result.add('cohort_disk_attention')
  if d['free_GiB'] < 60:result.add(d['host']+'_disk_attention')
 return result
state=r/'foreground-watch-seen.json';seen=set(json.loads(state.read_text())) if state.exists() else set()
start=time.monotonic();last_print=0
while True:
 with concurrent.futures.ThreadPoolExecutor(3) as ex:rows=list(ex.map(read,['aa-experiment','aa-screen','aa-nr']))
 found=events(rows)
 # Clear recovered alerts so a later recurrence wakes the foreground wait again.
 seen={key for key in seen if not any(tag in key for tag in ['job_error','connection_error','disk_attention','memory_attention','screen_failure']) or key in found}
 state.write_text(json.dumps(sorted(seen)))
 new=found-seen
 snapshot={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workers':rows,'new_events':sorted(new)}
 (r/'watch-status-latest.json').write_text(json.dumps(snapshot,indent=2))
 with (r/'progress-history.jsonl').open('a') as f:f.write(json.dumps(snapshot)+'\n')
 if new or time.monotonic()-last_print>300:
  keys=['host','native','repeat','primary_screens','screened','native_analysis_done','validation_started','validation_done','nr_reference_block','nr_reference_blocks','nr_done','nr_raw_backup_done','nr_memory_available_GiB','nr_swap_used_GiB','free_GiB']
  display={'utc':snapshot['utc'],'workers':[{k:d[k] for k in keys if k in d}|({'error_logs':d['error_logs']} if d.get('error_logs') else {})|({'connection_error':d['connection_error']} if d.get('connection_error') else {}) for d in rows if not d.get('expected_state')],'new_events':sorted(new)}
  print(json.dumps(display),flush=True);last_print=time.monotonic()
 if new:
  state.write_text(json.dumps(sorted(seen|new)));break
 if time.monotonic()-start>3600:break
 time.sleep(30)

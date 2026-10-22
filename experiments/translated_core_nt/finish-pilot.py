"""Foreground pilot sequence; wait for original run, compare frozen alternatives, then analyze."""
import datetime,json,subprocess,sys,time
from pathlib import Path
R=Path('/home/ec2-user/experiment');E=R/'experiments/translated_core_nt';D=R/'cohort-storage/translated-core-nt/pilot24'
commit=sys.argv[1]
plan=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_commit=commit,methods=['megablast','blastn','tblastx','dc-megablast','query-tblastn-fast','tblastx-word4'],policy='Sequential full-reference searches; identical original query panel and reference; original tblastx must complete before alternatives; preserve all raw evidence; no comparison on partial output')
(D/'alternative-sequence-protocol.json').write_text(json.dumps(plan,indent=2)+'\n')
while not (D/'tblastx/complete.json').exists():
 p=D/'tblastx/resources.json'
 if p.exists() and json.loads(p.read_text())['status']!='complete':raise RuntimeError('Original tblastx did not complete')
 if datetime.datetime.now(datetime.timezone.utc)>=datetime.datetime.fromisoformat('2026-10-23T08:00:00+00:00'):raise RuntimeError('Experiment cutoff')
 time.sleep(30)
for name in plan['methods'][3:]:
 if (D/name/'complete.json').exists():raise RuntimeError('Unexpected existing completed method: '+name)
 print(json.dumps(dict(starting_method=name,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())),flush=True)
 cmd=['python3',str(E/'run-query-only.py'),str(D),'--commit',commit] if name.startswith('query-') else ['python3',str(E/'experiment.py'),'run','--directory',str(D),'--config',str(E/(name+'.json')),'--commit',commit]
 subprocess.run(cmd,check=True)
subprocess.run(['python3',str(E/'analyze.py'),str(D),'--methods']+plan['methods'],check=True)
subprocess.run(['python3',str(E/'audit-codons.py'),'--directory',str(D)],check=True)
(D/'pilot-sequence-complete.json').write_text(json.dumps(dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_commit=commit,methods=plan['methods']),indent=2)+'\n')
print('PILOT_ALTERNATIVES_AND_ANALYSIS_COMPLETE',flush=True)

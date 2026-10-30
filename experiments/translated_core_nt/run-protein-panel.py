"""Foreground, bounded matched searches on the frozen larger diagnostic panel."""
import argparse,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--commit',required=True);a=p.parse_args();assert re.fullmatch('[a-f0-9]{40}',a.commit)
e=Path(__file__).resolve().parent;d=Path('/home/ec2-user/experiment/cohort-storage/translated-core-nt/protein232')
if not d.exists():subprocess.run(['python3',str(e/'prepare-protein-panel.py')],check=True)
plan=json.loads((e/'protein232-plan.json').read_text())
for method,seconds in zip(plan['methods_in_order'],plan['maximum_seconds']):
 print('START '+method,flush=True)
 subprocess.run(['python3',str(e/'experiment.py'),'run','--directory',str(d),'--config',str(e/(method+'.json')),'--commit',a.commit,'--seconds',str(seconds)],check=True)
 print('COMPLETE '+method,flush=True)
print('MATCHED_PROTEIN232_SEARCHES_COMPLETE',flush=True)

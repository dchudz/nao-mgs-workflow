#!/bin/bash
set -euo pipefail
cd /home/ec2-user/experiment
python3 - <<'PY'
import os,subprocess
from pathlib import Path
assert os.cpu_count()==32
assert int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemTotal:')))>=240*1024**2
assert subprocess.check_output(['aws','--profile','pipeline-codex','sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
assert Path('protein/nr24-oom-20261022T034211/attempt.json').exists()
assert not Path('protein/combined-validation/complete.json').exists()
PY
swapfile=/var/lib/aa-experiment-nr.swap
if ! sudo test -e "$swapfile"; then
 sudo fallocate -l 64G "$swapfile"
 sudo chmod 600 "$swapfile"
 sudo mkswap "$swapfile" >/dev/null
fi
if ! sudo swapon --show=NAME --noheadings | grep -Fxq "$swapfile"; then
 sudo swapon "$swapfile"
fi
sudo sysctl vm.swappiness=1
bash install-budget-safeguard.sh
python3 - <<'PY'
import datetime,json,os,subprocess
from pathlib import Path
info={line.split(':',1)[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemTotal:','MemAvailable:','SwapTotal:','SwapFree:'))}
assert info['SwapTotal']>=63*1024**2
record={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CPU_count':os.cpu_count(),'memory_KiB':info,'timers':subprocess.check_output(['systemctl','list-timers','--all','--no-pager'],text=True)}
assert 'Fri 2026-10-23 10:00:00 UTC' in record['timers']
Path('nr-memory-restart-host.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))
PY

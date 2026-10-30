"""Check real receipt acceptance and rejection of missing/stale evidence."""
import json
import shutil
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch
from synthetic_control_gate import verify_control

root=Path(sys.argv[1])
source=root/'source';control=root/'control-acceptance-v1'
diamond=Path('/home/ec2-user/experiment/protein/tools/diamond')
image='pipeline-e2e-rust:local'
proof=verify_control(source,diamond,image,control)
rejected=[]
def reject(name,call):
    try:call()
    except (ValueError,FileNotFoundError):rejected.append(name)
    else:raise AssertionError('Gate accepted '+name)
with tempfile.TemporaryDirectory(prefix='synthetic-gate-') as directory:
    tmp=Path(directory)
    reject('missing_receipt',lambda:verify_control(source,diamond,image,tmp))
    result=json.loads((control/'synthetic-control-results/attempt-3/result.json').read_text())
    for name in result['source_hashes']:
        dest=tmp/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source/name,dest)
    changed=tmp/'source/experiments/amino_acids/diamond_adapter.py'
    changed.write_text(changed.read_text()+'\n# simulated source change\n')
    reject('changed_adapter',lambda:verify_control(tmp/'source',diamond,image,control))
    shutil.copytree(control,tmp/'control')
    result['passed']=False
    (tmp/'control/synthetic-control-results/attempt-3/result.json').write_text(json.dumps(result))
    reject('failed_receipt',lambda:verify_control(source,diamond,image,tmp/'control'))
    with patch('synthetic_control_gate.subprocess.check_output',return_value='sha256:changed\n'):
        reject('changed_container',lambda:verify_control(source,diamond,image,control))
assert len(rejected)==4
record={'passed':True,'accepted_current_source_and_container':proof,'rejected':rejected}
(control/'synthetic-control-results/gate-verification.json').write_text(json.dumps(record,indent=2)+'\n')
proof['retrospective']=True
proof['note']='The first comparison was already launched before this control was requested.'
(root/'comparisons/SRR38294967_5m/synthetic-control-gate.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(record,indent=2))

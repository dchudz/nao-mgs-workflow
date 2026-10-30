"""Preview workflow wiring only in isolated output/work directories."""
import hashlib,json,subprocess
from pathlib import Path
R=Path('/home/ec2-user/experiment');D=R/'integration/parallel-adapter/preview';D.mkdir(exist_ok=True)
active=R/'e2e/SRR37076045/protein-run/trace.tsv';before=active.read_bytes()
params=json.loads((active.parent/'params.json').read_text());params['base_dir']=str(D)
(D/'params.json').write_text(json.dumps(params,indent=2)+'\n')
command=[str(R/'nextflow'),'-log',str(D/'nextflow.log'),'run',str(R),'-main-script','experiments/amino_acids/protein-parallel.nf','-lib',str(R/'lib'),'-c','experiments/amino_acids/e2e-worker.config','-c','e2e/SRR37076045/extra.config','-params-file',str(D/'params.json'),'-work-dir',str(D/'work'),'-preview','-ansi-log','false']
with (D/'preview.log').open('w') as f:subprocess.run(command,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
assert active.read_bytes().startswith(before),'Preview changed the active trace prefix'
assert len((D/'trace.tsv').read_text().splitlines())==1,'Preview executed a process'
receipt={'verified':True,'command':command,'active_trace_prefix_sha256':hashlib.sha256(before).hexdigest(),'active_trace_prefix_unchanged':True,'executed_processes':0,'output_directory':str(D)}
(D/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

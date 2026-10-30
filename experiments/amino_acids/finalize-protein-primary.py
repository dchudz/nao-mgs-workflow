"""Join completed comprehensive NR evidence to the frozen primary comparisons."""
import datetime, json, subprocess, time
from pathlib import Path
r=Path('/home/ec2-user/experiment'); c=r/'cohort-storage'; p=c/'protein200'; nr=r/'protein/combined-validation'; nr.mkdir(parents=True,exist_ok=True)
base='s3://pipeline-lab-839993661825-us-east-1/runs/amino-acid-20261021/'
aws=['aws','--profile','pipeline-codex','--region','us-east-1'];stop=datetime.datetime.fromisoformat('2026-10-23T09:00:00+00:00')
assert subprocess.check_output(aws+['sts','get-caller-identity','--query','Account','--output','text'],text=True).strip()=='839993661825'
logs=c/'protein-finalization'; logs.mkdir(exist_ok=True)
def run(script,*args):
 command=['python3',str(r/script),*map(str,args)]
 with (logs/'commands.jsonl').open('a') as f:f.write(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':command})+'\n')
 print('START '+script,flush=True)
 index=len((logs/'commands.jsonl').read_text().splitlines())
 with (logs/(str(index)+'-'+Path(script).name+'.log')).open('w') as f:subprocess.run(command,cwd=r,stdout=f,stderr=f,check=True)
 print('COMPLETE '+script,flush=True)
def backup(directory,prefix):
 subprocess.run(aws+['s3','cp',str(directory)+'/',base+prefix+'/', '--recursive','--only-show-errors'],check=True)
 run('verify-private-prefix.py',directory,'runs/amino-acid-20261021/'+prefix,'--output',directory/'backup-verification.json')
while True:
 result=subprocess.run(aws+['s3api','head-object','--bucket','pipeline-lab-839993661825-us-east-1','--key','runs/amino-acid-20261021/protein-nr-combined/complete.json'],capture_output=True,text=True)
 if result.returncode==0:break
 if '(404)' not in result.stderr and 'Not Found' not in result.stderr:raise RuntimeError(result.stderr)
 assert datetime.datetime.now(datetime.timezone.utc)<stop,'No complete NR evidence before the financial safeguard'
 time.sleep(60)
for name in ['complete.json','validated-pairs.tsv','native-control-validation.tsv','native-control-validation.json','seed-validation.tsv','seed-validation.json','search-protocol.json']:
 subprocess.run(aws+['s3','cp',base+'protein-nr-combined/'+name,str(nr/name),'--only-show-errors'],check=True)
protocol=json.loads((nr/'search-protocol.json').read_text());assert protocol['queries']==817744 and protocol['query_sha256']=='6aa28bedbaf612c3a7c5a895b24a3be4f18d2c8264a439846427d943b289c597'
for arm,suffix in [('nt',''),('nt20','-nt20')]:
 out=p/('native-comparison'+suffix)
 run('analyze-protein200.py','--protein',p,'--nr',nr/'validated-pairs.tsv','--native-arm',arm,'--output',out)
 run('analyze-protein-overlap.py',out/'pair-evidence.tsv',out/'original-read-overlap.json','--gate',r/'protein-repeat-gate/query-gate.tsv','--screen-membership',r/'protein-screen-boundaries/query-membership.tsv','--native-arm',arm)
 assert json.loads((out/'summary.json').read_text())['nr_complete']
 overlap=json.loads((out/'original-read-overlap.json').read_text());assert 'pending' not in overlap['groups'] and overlap['repeat_gate_included'] and overlap['initial_screen_included']
 backup(out,'protein200-native-comparison'+suffix)
run('analyze-native-protein-losses.py','--nr-controls',nr/'native-control-validation.tsv','--raw-screen',r/'protein-repeat-gate/native-loss-raw-screen.tsv')
backup(p/'native-losses','protein-native-losses')
run('experiments/amino_acids/compare-protein.py','--metadata',r/'full/protein/queries.tsv.gz','--members',r/'full/protein/queries.fna.members.tsv.gz','--viral',r/'full/protein/viral-pairs.tsv','--nr',nr/'validated-pairs.tsv','--native',r/'full/control/postscreen/RENAME_VIRUS_HITS/SRR39694979_virus_hits.tsv.gz','--fastp',r/'full/control/postscreen/FASTP/SRR39694979_fastp.fastq.gz','--prefix','FULL_','--output',r/'full/protein/native-comparison.tsv')
full=logs/'full-anchor';full.mkdir(exist_ok=True)
for name in ['native-comparison.tsv','native-comparison.tsv.json']:(full/name).write_bytes((r/'full/protein'/name).read_bytes())
backup(full,'full-protein-native-comparison')
run('summarize-protein-nt-panel.py',c/'protein-nt-panel','--nodes',r/'reference/results/taxonomy-nodes.dmp','--nr',nr/'validated-pairs.tsv')
assert json.loads((c/'protein-nt-panel/comparison.json').read_text())['nr_complete']
backup(c/'protein-nt-panel','protein-nt-panel')
run('compare-validation-methods.py','--nucleotide',c/'seed-diagnostic200/comparison.tsv','--protein',nr/'seed-validation.tsv','--output',c/'validation-methods')
backup(c/'validation-methods','validation-methods')
marker={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protein_samples':200,'native_baselines':['nt','nt20'],'full_depth_anchor':True,'panel_pairs':232,'diagnostic_pairs':746,'private_backups_verified':True}
(logs/'complete.json').write_text(json.dumps(marker,indent=2))
subprocess.run(aws+['s3','cp',str(logs)+'/',base+'protein-finalization/','--recursive','--exclude','full-anchor/*','--exclude','complete.json','--only-show-errors'],check=True)
subprocess.run(aws+['s3','cp',str(logs/'complete.json'),base+'protein-finalization/complete.json','--only-show-errors'],check=True)
print('PRIMARY_PROTEIN_ANALYSES_COMPLETE',flush=True)

"""Back up local experiment history and evidence using in-memory worker credentials."""
import datetime,hashlib,json,os,pathlib,subprocess,tarfile,tempfile
root=pathlib.Path('/Users/david/Documents/ChatGPT/pipeline');base=root/'.tmp';d=base/'translated-core-nt-20261022';work=pathlib.Path('/Users/david/.codex/worktrees/amino-acid-screen/pipeline')
head=subprocess.check_output(['git','-C',str(work),'rev-parse','HEAD'],text=True).strip();short=head[:12];old=json.loads((d/'latest-backup-manifest.json').read_text())
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
retirement_args=['--retirement-receipts',str(d/'owned-worker-retirements.json')] if (d/'owned-worker-retirements.json').exists() else []
subprocess.run(['python3',str(base/'aa-experiment-20261021/estimate-experiment-cost.py'),'--as-of',now,'--output',str(d/'budget-ledger.json')]+retirement_args,check=True,stdout=subprocess.DEVNULL)
ledger=json.loads((d/'budget-ledger.json').read_text());split=json.loads((d/'split-budget-ledger.json').read_text());split['as_of_utc']=now
for k,v in split['before_boundary'].items():split['after_boundary'][k]=ledger[k]-v
split['new_cap_modeled_remainder_excluding_other_charges']=500-split['after_boundary']['modeled_owned_compute_EBS_IPv4_total_USD'];(d/'split-budget-ledger.json').write_text(json.dumps(split,indent=2)+'\n')
state=json.loads((d/'active-runtime.json').read_text());state.update(updated_utc=now,latest_budget=split,current_local_head=head,backup_refresh_needed=True);(d/'active-runtime.json').write_text(json.dumps(state,indent=2)+'\n')
bundle=d/f'history-{short}.bundle';source=d/f'source-{short}.tar.gz';reports=d/f'reports-{short}.tar.gz';assert not any(p.exists() for p in [bundle,source,reports])
subprocess.run(['git','-C',str(work),'bundle','create',str(bundle),'aa3dcb073af33f6331f6c4223c656e4bab7cd2eb..HEAD'],check=True)
subprocess.run(['git','-C',str(work),'bundle','verify',str(bundle)],check=True,stdout=subprocess.DEVNULL)
subprocess.run(['git','-C',str(work),'archive','--format=tar.gz','--output='+str(source),head,'experiments/amino_acids','experiments/translated_core_nt'],check=True)
prior=next(x['name'] for x in old['files'] if x['name'].startswith('reports-'))
with tarfile.open(d/prior,'r:gz') as t:members={x.name for x in t.getmembers() if x.isfile()}
snapshot=d/f'foreground-progress-{short}.jsonl';snapshot.write_bytes((d/('foreground-replay-progress-20261023.jsonl' if (d/'foreground-replay-progress-20261023.jsonl').exists() else 'foreground-progress-20261023.jsonl')).read_bytes())
add=[d/'protein232',d/'protein232-result-v2.html',d/'protein232-result-v2.provenance.json',d/'protein232-download-verification.json',d/'measured-search-costs-v3',d/'watch-final-comparisons.py',d/'blast-2.17-batching-source',d/'measured-search-costs-v2',d/'protein232/query-batching-release-verification.json',d/'protein232/query-batching-progress-observation.json',d/'pulls-912-public-current.json',d/'issues-920-public-current.json',d/'budget-ledger.json',d/'split-budget-ledger.json',d/'active-runtime.json',snapshot,d/'batch-replay-characterization.json',d/'watch-active-comparisons.py',pathlib.Path(__file__)]
add += [d/'seed-nr-forced-double-index-trace',d/'fixed-nr-audit-download-verification.json',d/'nr-context-result-v2.html',d/'nr-context-result-v2.provenance.json',d/'query-fast-stop-score-characterization.json',d/'measured-search-costs-v4',d/'nr-fixed-mode-ssh-transport-disconnection.json']
add += [base/'aa-experiment-20261021/nr200-result-v3.html',base/'aa-experiment-20261021/nr200-result-v3.provenance.json',d/'nr-worker-retirement.json',d/'owned-worker-retirements.json',d/'nr-final-evidence-complete.json',d/'nr-reference-archive-complete.json',d/'nr-reference-preservation-summary.json',d/'retire-preserved-nr-worker.py',d/'measured-search-costs-v5',d/'pair-analysis-interface-parity.json',d/'protein232-original-panel-current',d/'dc-megablast-protein232-plan.json']
add += [d/'canonical-score-disagreement-diagnostic.json', d/'diamond-2.1.17-source/tree/src/stats/matrices/blosum62.h', d/'one-batch-final-receipt.json', d/'coding-check-not-run.json', d/'protein232-result-v3.html', d/'protein232-result-v3.provenance.json', d/'protein232-report-v3-verification.json', d/'measured-search-costs-v6', d/'quota-verification-20261023.json', d/'preservation-safeguard-final-applied.json', d/'cost-explorer-latest-summary.json']
for optional in ['primary-final-evidence', 'primary-worker-retirement.json', 'core-reference-archive-complete.json', 'retained-storage-final.json', 'experiment-results.html', 'experiment-results.provenance.json', 'write-experiment-results.py', 'retire-preserved-primary-worker.py', 'primary-symlink-coverage.json', 'refresh-billing-storage.py']:
 if (d/optional).exists():add.append(d/optional)
for p in add:
 assert p.exists(),p
 for f in p.rglob('*') if p.is_dir() else [p]:
  if f.is_file():members.add(str(f.relative_to(base)))
with tarfile.open(reports,'w:gz') as t:
 for name in sorted(members):
  p=base/name;assert p.is_file();t.add(p,arcname=name,recursive=False)
def record(p):return {'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
m={'head':head,'created_utc':now,'report_files':len(members),'files':[record(p) for p in [bundle,source,reports]]}
(d/'latest-backup-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
# Credentials are captured in memory only and never written into a file or output.
if os.environ.get('PIPELINE_CLOSURE_USE_INHERITED_AWS_CREDENTIALS') == '1':
 assert all(os.environ.get(k) for k in ['AWS_ACCESS_KEY_ID','AWS_SECRET_ACCESS_KEY','AWS_SESSION_TOKEN'])
 env=os.environ.copy()
else:
 creds=json.loads(subprocess.check_output(['ssh','-F',str(base/'aa-experiment-20261021/ssh_config'),'aa-experiment','aws configure export-credentials --profile pipeline-codex --format process'],text=True))
 env={k:v for k,v in os.environ.items() if not k.startswith('AWS_')};env.update(AWS_ACCESS_KEY_ID=creds['AccessKeyId'],AWS_SECRET_ACCESS_KEY=creds['SecretAccessKey'],AWS_SESSION_TOKEN=creds['SessionToken'])
env.update(AWS_CONFIG_FILE='/dev/null',AWS_SHARED_CREDENTIALS_FILE='/dev/null',AWS_EC2_METADATA_DISABLED='true',AWS_DEFAULT_REGION='us-east-1')
a=['/opt/homebrew/bin/aws','--region','us-east-1'];assert subprocess.check_output(a+['sts','get-caller-identity','--query','Account','--output','text'],env=env,text=True).strip()=='839993661825'
prefix='s3://pipeline-lab-839993661825-us-east-1/runs/translated-core-nt-20261022/source-revisions/'+head+'/'
def upload_verify(p):
 expected=record(p);subprocess.run(a+['s3','cp',str(p),prefix+p.name,'--only-show-errors'],env=env,check=True)
 with tempfile.TemporaryDirectory() as tmp:
  q=pathlib.Path(tmp)/p.name;subprocess.run(a+['s3','cp',prefix+p.name,str(q),'--only-show-errors'],env=env,check=True);assert record(q)==expected
for p in [bundle,source,reports]:upload_verify(p)
receipt={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'account':'839993661825','prefix':prefix,'head':head,'files':m['files'],'all_download_hashes_verified':True}
(d/'latest-source-backup-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
state.update(last_verified_source_backup=m,backup_refresh_needed=False);(d/'active-runtime.json').write_text(json.dumps(state,indent=2)+'\n')
for p in [d/'latest-backup-manifest.json',d/'active-runtime.json',d/'latest-source-backup-verification.json']:upload_verify(p)
print(json.dumps(receipt),flush=True)

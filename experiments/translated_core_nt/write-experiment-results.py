"""Link the finished private experiments and their preservation/cost receipts."""
import datetime,hashlib,html,json,os
from pathlib import Path
p=Path(__file__).resolve().parent
old=p.parent/'aa-experiment-20261021'
inputs={}
def read(path):
 data=path.read_bytes();inputs[str(path)]=hashlib.sha256(data).hexdigest();return json.loads(data)
def link(path,label):
 assert path.exists(),path
 return '<a href="'+html.escape(os.path.relpath(path,p))+'">'+html.escape(label)+'</a>'
split=read(p/'split-budget-ledger.json')
quota=read(p/'quota-verification-20261023.json')
batch=read(p/'one-batch-final-receipt.json')
state=read(p/'active-runtime.json')
retired=read(p/'owned-worker-retirements.json')
rows=[('First upstream step','1,000 samples; first five million pairs per sample',old/'first-step-result.html'),('Protein detection and NR checking','200-sample subset; experimental call generation',old/'nr200-result-v3.html'),('Translated search against the same core_nt','232 selected pairs / 464 mates',p/'protein232-result-v3.html'),('NR search-mode reproducibility','746 selected pairs / 1,436 mates',p/'nr-context-result-v2.html'),('Translated search algorithms','43-mate pilot and scoring inspection',p/'faster-translated-search-v3.html'),('Read-header reproducibility','Eight selected read IDs',old/'bowtie-header-note.html')]
text='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Amino-acid experiment results</title><style>body{font:17px/1.5 system-ui;max-width:1080px;margin:40px auto;padding:0 24px;color:#172332}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #cbd3dc;padding:10px;text-align:left}a{color:#165d9c}.note{background:#edf2f7;padding:16px}code{overflow-wrap:anywhere}</style><h1>Amino-acid experiment results</h1><p>Private results. AI-authored analysis, not human-verified. No source or results were pushed or publicly posted.</p>'''
text+='<p class="note">These are completed component experiments and explicitly bounded diagnostic attempts. A full production amino-acid pipeline over the entire public dataset was not run. Each linked report states its input selection, filters, measured results and limitations.</p>'
text+='<table><tr><th>Experiment</th><th>Scope</th><th>Evidence</th></tr>'
for title,scope,path in rows:text+='<tr><td>'+title+'</td><td>'+scope+'</td><td>'+link(path,'Report and exact data links')+'</td></tr>'
score=read(p/'canonical-score-disagreement-diagnostic.json')
assert score['counts']['canonical/different']==6219 and score['counts']['canonical/equal']==4784
text+='</table><p>An additional score audit found different bit scores for 6,219 of 11,003 identical shared alignments containing only the 20 standard amino acids. Removing stop characters alone therefore does not establish score equivalence. The retained output does not contain enough scoring-context data to reconstruct or calibrate those scores. '+link(p/'canonical-score-disagreement-diagnostic.json','Canonical-amino-acid scoring diagnostic')+'.</p>'
text+='<h2>Final bounded attempt</h2><p>The one-batch translated replay ended at its preset time limit after '+f"{batch['search_status']['wall_seconds']:,.1f}"+' seconds, with no finalized output. Its modeled search allocation was $'+f"{batch['allocated_compute_EBS_IPv4_USD']:.2f}"+'. It does not establish speedup, output parity or biological no-hits. The optional coding-template check was not launched because its planned launch window had passed. '+link(p/'one-batch-final-receipt.json','Replay receipt')+'; '+link(p/'coding-check-not-run.json','Unexecuted plan receipt')+'.</p>'
text+='<h2>Spend and retained evidence</h2><p>Modeled EC2, EBS and IPv4 allocation: $'+f"{split['before_boundary']['modeled_owned_compute_EBS_IPv4_total_USD']:.2f}"+' before the October 22 14:51 UTC budget boundary; $'+f"{split['after_boundary']['modeled_owned_compute_EBS_IPv4_total_USD']:.2f}"+' after it. Each window has its own $500 cap. This model includes staging, interrupted searches, analysis, archival and idle time. Per-search costs are already included. S3, requests, transfer, tax and credits are separate. '+link(p/'split-budget-ledger.json','Separate budget ledgers')+'.</p>'
if (p/'retained-storage-final.json').exists():
 storage=read(p/'retained-storage-final.json')
 text+='<p>Retained private S3 objects are modeled at $'+f"{storage['modeled_monthly_data_storage_USD']:.2f}"+'/month. Archived references require restoration before reuse. Request, retrieval, transfer, version and incomplete-upload charges are excluded. '+link(p/'retained-storage-final.json','Storage inventory and rates')+'.</p>'
text+='<p>Worker retirement receipts: '+link(p/'nr-worker-retirement.json','NR worker')
if (p/'primary-worker-retirement.json').exists():text+='; '+link(p/'primary-worker-retirement.json','Primary worker')
else:text+='; primary preservation is still in progress'
text+='.</p>'
for name,label in [('core-reference-archive-complete.json','Exact core_nt archive'),('nr-reference-archive-complete.json','Exact NR archive'),('primary-final-evidence/complete.json','Primary scientific archives'),('nr-final-evidence-complete.json','NR scientific archive')]:
 if (p/name).exists():text+='<p>'+link(p/name,label)+'</p>'
text+='<p>The live EC2 API confirmed '+str(int(quota['quota']['Value']))+' standard on-demand vCPUs in us-east-1. No additional compute was launched following that quota increase. '+link(p/'quota-verification-20261023.json','Quota receipt')+'.</p>'
text+='<p>Experiment source, failed attempts, results and analysis revisions remain in separate local commits on <code>codex/amino-acid-screen</code>. '+link(p/'latest-source-backup-verification.json','Latest verified private history/source/report backup')+'.</p></html>'
out=p/'experiment-results.html';out.write_text(text)
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':inputs,'report_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'retired_instances':retired}
out.with_suffix('.provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(str(out))

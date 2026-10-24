"""Render saved four-arm results without changing any calling or selection rule."""
import argparse,csv,html,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--plan',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
esc=lambda x:html.escape(str(x))
labels={'production':'Production','validation_control':'Sensitive validation','screen_control':'Expanded screen + DNA alignment','protein_variant':'Expanded screen + protein alignment'}
host={'0':'Inconsistent','1':'Match','2':'Unclear','3':'Consistent','-1':'Unresolved','-2':'Maybe inconsistent','NA':'Unavailable','':'Unavailable'}
selection=json.loads((a.plan/'sample-selection.json').read_text())['initial_samples']
budget=json.loads((a.plan/'budget.json').read_text())
summary_path=a.root/'analysis/summary.json'
s=json.loads(summary_path.read_text()) if summary_path.exists() else {'completed_samples':[],'metrics':[],'changes':[],'pair_contrasts':[]}
completed=set(s['completed_samples'])
def table(headers,rows):
 rows=list(rows)
 if not rows:return '<p>Pending: no completed comparison.</p>'
 return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def integer(value):return f'{int(value):,}'
metrics=table(['Sample','Arm','Candidate pairs','Joint duplicate exemplars','Candidate species calls','Validated species calls'],[[x['sample'],labels[x['arm']],integer(x['rows']),integer(x.get('joint_duplicate_exemplars',0)),x['candidate_species_calls'],x['validated_species_calls']] for x in s['metrics']])
changes=table(['Sample','Variant / baseline','Tier','Change','Species (taxid)','Human-host annotation'],[[x['sample'],labels[x['variant']]+' / '+labels[x['baseline']],x['tier'],x['change'],x['name']+' ('+str(x['species_taxid'])+')',host.get(x['human_annotation'],x['human_annotation'])] for x in s['changes']]) if s['changes'] else '<p>'+('No species-call changes in the completed comparisons.' if completed else 'Pending: no completed comparison.')+'</p>'
pairs=table(['Sample','Variant / baseline','Candidate pairs gained','Lost','Shared'],[[x['sample'],labels[x['variant']]+' / '+labels[x['baseline']],integer(x['candidate_pairs_gained']),integer(x['candidate_pairs_lost']),integer(x['candidate_pairs_shared'])] for x in s['pair_contrasts']])
cohort=table(['Sample','Full-depth input pairs','Status'],[[x['run_accession'],integer(x['read_count']),'All four outputs complete' if x['run_accession'] in completed else 'Pending'] for x in selection])
timings=[];versions=[]
for row in selection:
 sample=row['run_accession'];d=a.root/sample
 if not d.exists():continue
 for stage in ['baseline-run','aa-screen','union','screen-run','protein-run','production-downstream','validation_control-downstream','screen_control-downstream','protein_variant-downstream']:
  f=d/(stage+'.json')
  if f.exists():
   r=json.loads(f.read_text());timings.append([sample,stage,f"{r['wall_seconds']/60:.2f}",r.get('exit_code','')]);versions.append((r.get('source_commit',''),r.get('driver_commit','')))
version_html=''.join('<li>Scientific source <code>'+esc(src)+'</code>; execution driver <code>'+esc(driver)+'</code></li>' for src,driver in sorted(set(versions)))
status_path=a.root/'worker-status.json'
status=json.loads(status_path.read_text()) if status_path.exists() else None
cost=(f"Worker allocation upper bound: ${status['worker_cost_upper_USD']:.2f} at {esc(status['utc'])}. This excludes transfer, retained object storage and requests; it is not the total AWS bill." if status else 'Worker-cost reconciliation pending.')
body=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Full-depth pipeline comparison</title><style>body{{font:16px/1.55 system-ui,sans-serif;color:#17212b;max-width:1400px;margin:40px auto;padding:0 24px}}h1,h2{{line-height:1.2}}h2{{margin-top:34px}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{border:1px solid #d4dde5;padding:8px;text-align:left;vertical-align:top}}th{{background:#edf2f6}}.scroll{{overflow-x:auto}}code{{overflow-wrap:anywhere}}.note{{background:#f2f5f8;padding:16px}}</style></head><body>
<h1>Full-depth pipeline comparison</h1><p>{len(completed)} of {len(selection)} frozen samples have completed all four outputs. Planned input: {integer(sum(int(x['read_count']) for x in selection))} read pairs. Partial stages appear only in the runtime table; they do not enter biological comparisons.</p>
<h2>Four outputs from the same full-depth inputs</h2><ol><li><strong>Production:</strong> NT24 screen, Bowtie2 alignment, default megablast validation.</li><li><strong>Sensitive validation:</strong> the same native candidates, validated with blastn -task blastn.</li><li><strong>Expanded screen + DNA alignment:</strong> NT24 OR AA12 screen, Bowtie2 alignment, sensitive blastn validation.</li><li><strong>Expanded screen + protein alignment:</strong> the same union and FASTP input as output 3, DIAMOND blastx alignment, sensitive blastn validation.</li></ol>
<p>The alignment reference retains the same DUST-masked viral genomes; the protein index translates them in six frames. Every validation arm uses the same comprehensive core_nt. The AA screen uses a human-masked index, so the screening contrast includes a masking difference. FASTP, human/other filtering, production duplicate marking and per-species/group validation cap20 are retained. Each sample forms its own group. QC and profiling are reused across arms because they consume the same raw input independently of viral screening; profiling uses the production one-million-pair subsample.</p>
<h2>Completed sample results</h2>{metrics}<p>A candidate species call requires at least three joint alignment-and-similarity duplicate exemplars. A validated call requires at least three such exemplars, each with its own aligned validation natural LCA resolving to the same species as its viral alignment. Unsampled pairs receive no inferred validation verdict. These are computational evidence criteria, not biological ground truth or evidence of infection.</p>
<h2>Species calls gained and lost</h2>{changes}<p>Human-host labels are copied from the frozen reference's infection_status_human codes. They do not establish infection in a sample. Artificial-reference competition is retained in the read-level evidence and separately counted in species-support.tsv; the frozen natural-LCA calling rule is unchanged.</p>
<h2>Candidate-pair differences</h2>{pairs}
<h2>Runtime and budget</h2><p>Fresh authorized cap: ${budget['cap_USD']}; no previous budget remainder carries forward. {cost}</p>{table(['Sample','Completed stage','Wall minutes','Exit code'],timings)}<p>Stage wall times include contention and differing cache state. Some stages overlap, so adding their wall times does not give billed worker time. Database staging, acquisition, preservation and analysis are additional work.</p>
<h2>Frozen cohort</h2>{cohort}
<h2>Evidence and provenance</h2><ul>{version_html}</ul><p>Analysis source: <code>{esc(s.get('analysis_source_commit','Pending'))}</code>. Full production merged/joined validation queries are used; this differs from the earlier NR both-mate criterion.</p><ul><li><a href="analysis/summary.json">Summary, definitions and limits</a></li><li><a href="analysis/species-support.tsv">All species support, including artificial-reference evidence</a></li><li><a href="analysis/species-changes.tsv">Species-call differences</a></li><li><a href="analysis/pair-changes.tsv">Candidate-pair differences</a></li><li><a href="analysis/read-disagreements.tsv.gz">Read-level disagreements</a></li><li><a href="analysis/call-evidence.tsv.gz">Individual exemplar evidence</a></li></ul><p>Complete commands, resource traces, per-read validation outputs, HSPs and diagnostic archives are retained in private sample checkpoints. Each completed sample must pass independent S3 readback verification before scratch is released. No results are published externally.</p></body></html>'''
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(body);print(json.dumps({'report':str(a.output),'completed_samples':len(completed)}))

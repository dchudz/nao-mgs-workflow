"""Render completed NR context measurements and clearly distinguish the pending mode control."""
import argparse,csv,html,json,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();r=a.root
load=lambda p:json.loads(p.read_text())
context=load(r/'seed-nr-standalone/context-cost-and-mode.json');c=context['reproducibility'];counts=load(r/'NR-context-competitor-analysis-v2/summary.json')['counts']
esc=lambda x:html.escape(str(x))
def table(headers,rows):
    return '<table><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'
def link(path,label):return '<a href="'+esc(os.path.relpath(path,a.output.parent))+'">'+esc(label)+'</a>'
search_rows=[]
for name,record in context['searches'].items():
    search_rows.append([name,','.join(record['algorithms']),f'{record["queries"]:,}',f'{record["GNU_wall_seconds"]:,.0f}',f'${record["allocated_compute_EBS_IPv4_USD"]:.3f}',f'${record["allocated_USD_per_query"]:.6f}'])
text='''<!doctype html><html lang="en"><meta charset="utf-8"><title>NR batch reproducibility</title><style>body{font:17px/1.55 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#172332}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #cbd3dc;padding:9px;text-align:left}a{color:#165d9c}.note{background:#edf2f7;padding:16px}code{overflow-wrap:anywhere}</style><h1>NR batch reproducibility</h1><p>Private experiment. AI-authored analysis, not human-verified.</p>'''
text+=f'<p><strong>{c["status_disagreements"]} of {c["pairs"]} diagnostic pairs change their original viral/nonviral status when the same query sequences are searched alone instead of inside the large pooled batch.</strong> Taxonomic assignments change for {c["taxid_disagreements"]} pairs. After distinguishing artificial constructs from natural competitors, {counts["natural_status_changes"]} pairs change natural-evidence status.</p>'
text+='<p class="note">This is a selected reproducibility diagnostic, separate from the full 200-sample candidate set. It has no known-truth labels. Neither disagreement fractions nor apparent improvements can be extrapolated to the candidate calls or a public population.</p>'
text+='<h2>Inputs and observed execution modes</h2><p>The same 1,436 mate sequences from 746 frozen pairs appear in both runs, with unchanged DIAMOND 2.1.17, comprehensive NR snapshot and search options. Both use sensitive mode, 32 threads, block size 2, four index chunks, minimum ORF 1, genetic code 1, E≤10⁻⁵, 30% query coverage and the top 10% score range. The postfilter requires at least 25 aligned amino acids and taxa within 90% of each mate’s best eligible score; pair taxa combine only that pair’s own mates.</p><p>DIAMOND automatically selected query-indexed search for the small input and double-indexed search for the large input. The <a href="https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options">official options documentation</a> describes this automatic choice and warns that sensitivity can differ. These completed runs change both query context and the automatically selected mode; they do not isolate those effects.</p>'
text+=table(['Context','Logged mode','Query mates','Search seconds','Allocated search cost','Cost/query'],search_rows)
text+='<p>Search allocation includes the 32-vCPU r7a.8xlarge, provisioned gp3 and public IPv4 over GNU elapsed search time, including the recorded volume expansion. It excludes analysis, backups, staging, idle time and S3. The lifetime ledger already includes these costs; do not add them again. Different query pools and cache states prevent an isolated algorithm-speed interpretation.</p>'
text+='<h2>Returned alignments and classifications</h2>'
text+=table(['Measurement','Count'],[['Pooled diagnostic alignment rows',f'{c["pooled_alignment_rows"]:,}'],['Standalone alignment rows',f'{c["standalone_alignment_rows"]:,}'],['Rows only in pooled output',f'{c["alignment_rows_only_pooled"]:,}'],['Rows only in standalone output',f'{c["alignment_rows_only_standalone"]:,}'],['Mates with raw-output differences',c['raw_alignment_differing_query_mates']],['Pairs with any summary-field difference',c['pairs_with_summary_difference']]])
text+=table(['Original status transition','Pairs'],counts['original_status_transitions'].items())
text+=f'<p>Of {counts["mixed_to_viral_pairs"]} pairs changing from mixed to viral-only, {counts["mixed_to_viral_with_prior_artificial_only_nonviral"]} previously had artificial constructs as their only nonviral competitors. Another {counts["mixed_to_viral_with_prior_cellular"]} had cellular competitors, including {counts["mixed_to_viral_with_prior_human"]} with human matches. Artificial constructs are not described as human or bacterial competition.</p>'
text+=table(['Natural-evidence transition','Pairs'],counts['changed_natural_status_transitions'].items())
text+='<h2>Two human-competition cases</h2><p>The human matches for <code>SRR33038061.678928042</code> and <code>SRR35987653.689121697</code> are absent from the standalone raw output, rather than being removed only by our downstream taxonomic filter. In both cases the original human scores would still exceed 90% of the standalone best score. A human match does not establish a human source for the read.</p>'
human=list(csv.DictReader((r/'NR-context-human-cases/human-HSPs.tsv').open(),delimiter='\t'))
text+=table(['Query','Human bit score','E-value','Aligned amino acids','Identity','Target label'],[[x['query'],x['bitscore'],x['evalue'],x['alignment_AA'],x['identity_pct']+'%',x['subject']] for x in human])
text+='<p>These four human HSPs have an unavailable accession label (<code>N/A</code>) in the retained original output. A separate synthetic native-BLAST fixture verified that adding the numeric target ordinal leaves all seven original alignment columns unchanged in both modes and allows exact target-sequence retrieval even when a label is missing. The control below retains that additional identifier; no target accession is inferred from a taxonomic label.</p>'
text+='<h2>Fixed-mode control</h2>'
fixed=r/'seed-nr-forced-double-index-trace'
if (fixed/'complete.json').exists():
    comparisons=load(fixed/'context-comparison.json')
    text+=table(['Comparison to fixed-mode standalone','Raw multiset equal','Status differences','Taxonomic differences','Any summary difference'],[[name,x['alignment_multiset_equal'],x['status_disagreements'],x['taxid_disagreements'],x['pairs_with_summary_difference']] for name,x in comparisons.items()])
    text+='<p>'+link(fixed/'context-comparison.json','Completed controlled comparison')+'</p>'
else:
    text+='<p>One repeat of the identical 1,436-query standalone set is running with double-indexed mode forced (<code>--algo 0</code>). It retains the original reference, executable, threads and search thresholds. The additional <code>snum</code> column records target provenance; a lossless seven-column projection feeds the original classifier. The run is capped at seven hours, at most $18.76 of allocated search cost. No result from this control is yet included above.</p>'
text+='<p>The first incomplete control attempt was stopped after 1,512.65 seconds ($1.13 allocated) to add target traceability; it produced no final raw alignment rows. Its protocol, log and interruption receipt are retained. It is excluded from completed-search comparisons and included in lifetime spend.</p>'
text+='<h2>Reference work and verification</h2><p>Both searches use the existing 918.57 GB native BLAST-format NR snapshot: 1,172,457,802 sequences and 441,882,733,350 amino acids. No persistent DIAMOND index was built. Reference seed structures are rebuilt during each run.</p>'
phase_rows=[]
for key in context['searches']['combined-validation']['logged_stage_groups']:
    row=[key.replace('_',' ')]
    for record in context['searches'].values():
        group=record['logged_stage_groups'][key];row.append(f'{group["logged_seconds"]:,.1f}' if group['labels'] else 'not separately logged')
    phase_rows.append(row)
text+=table(['Logged phase','Pooled elapsed seconds','Standalone elapsed seconds'],phase_rows)
text+='<p>These are logged elapsed stage durations, not CPU time; they need not sum to the overall timer. A stage not separately logged is not proof that no such work occurred. Local evidence annotation exactly reproduced all 1,492 saved context-specific pair classifications, taxon counts and best scores before adding competitor categories. Raw downloads were verified against source SHA-256 hashes.</p><h2>Evidence</h2>'
for path,label in [(r/'seed-nr-standalone/context-cost-and-mode.json','Commands, timing provenance and allocation'),(r/'seed-nr-standalone/context-disagreements.tsv','Every original pair-summary disagreement'),(r/'NR-context-competitor-analysis-v2/per-pair-competitors.tsv','Every pair’s artificial and cellular competitors'),(r/'NR-context-human-cases/query-comparison.tsv','Both raw-output human-competition cases'),(r/'NR-context-local-evidence/download-receipt.json','Raw-input hash verification'),(r/'nr-target-ordinal-fixture-v3/complete.json','Synthetic output and retrieval fixture'),(r.parent/'aa-experiment-20261021/nr200-result-v2.html','Original pooled candidate results with context limitation')]:
    assert path.exists(),path;text+='<p>'+link(path,label)+'</p>'
text+='</html>';a.output.write_text(text);print(a.output)

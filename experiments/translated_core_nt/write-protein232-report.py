"""Render completed same-reference evidence without assigning biological truth."""
import argparse
import collections
import csv
import datetime
import hashlib
import html
import json
import os
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('directory', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--costs', type=Path, required=True)
p.add_argument('--codon-audit', type=Path, required=True)
a = p.parse_args()
d = a.directory
methods = ['megablast', 'blastn', 'tblastx-word4']
inputs = {}
def read(path):
    data = path.read_bytes()
    inputs[str(path)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)
manifest = read(d / 'input-manifest.json')
assert manifest['pairs'] == 232 and manifest['queries'] == 464
for method in methods:
    assert read(d / method / 'resources.json')['status'] == 'complete'
    assert read(d / method / 'complete.json')['queries'] == 464
mate = read(d / 'analysis/summary.json')
tax = read(d / 'analysis/taxonomy-description/summary.json')
pair = read(d / 'pair-analysis/summary.json')
audit = read(a.codon_audit)
transitions = read(d / 'transition-analysis/summary.json')
routing = read(d / 'selective-validation/summary.json')
cost = read(a.costs)
assert mate['pairs'] == pair['pairs'] == 232 and mate['mates'] == pair['mates'] == 464
assert pair['NR_frozen_pair_verdict_parity_verified']
searches = {x['method']: x for x in cost['searches'] if x['panel'] == 'protein232' and x['method'] in methods}
assert set(searches) == set(methods)
assert all(x['query_sha256'] == manifest['query_sha256'] for x in searches.values())
esc = lambda x: html.escape(str(x))
def table(headers, rows):
    return '<table><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'
def link(path, label):
    assert path.exists(), path
    return '<a href="'+esc(os.path.relpath(path, a.output.parent))+'">'+esc(label)+'</a>'
def counts_table(groups, key):
    columns = sorted({k for group in groups.values() for k in group[key]})
    return table(['Group', 'Pairs']+columns, [[name, x['pairs']]+[x[key].get(k, 0) for k in columns] for name, x in sorted(groups.items())])

text='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Same-reference translated validation</title><style>body{font:17px/1.5 system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#172332}table{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0}td,th{border:1px solid #cbd3dc;padding:8px;text-align:left;overflow-wrap:anywhere}a{color:#165d9c}.note{background:#edf2f7;padding:16px}code{overflow-wrap:anywhere}details{margin:20px 0}</style><h1>Same-reference translated validation</h1><p>Private experiment. AI-authored analysis, not human-verified.</p>'''
text+='<p>Three completed searches compare the same 464 mate sequences from 232 frozen pairs against the pipeline’s comprehensive <code>core_nt</code> snapshot: megablast, sensitive nucleotide BLAST, and native translated BLAST with four-amino-acid seeds. Translation changes the search; it does not substitute a different reference database.</p>'
text+='<p class="note">The panel contains one protein-only candidate from each of 200 samples plus 32 shared controls. Selection was frozen before NR outcomes. This is a selected validation comparison, not a full pipeline run, a population accuracy estimate, or a demonstration of new sample/species calls at the three-exemplar threshold. Every read is judged from its own alignments.</p>'
text+='<h2>Measured search work</h2>'
text+=table(['Method','Query mates','Wall seconds','CPU seconds','Allocated USD','USD/query','Raw HSP rows'],[[m,searches[m]['query_mates'],f"{searches[m]['wall_seconds']:,.2f}",f"{searches[m]['CPU_seconds']:,.2f}",f"{searches[m]['allocated_USD']:.4f}",f"{searches[m]['allocated_USD_per_query']:.6f}",mate['resources'][m]['raw_alignment_rows']] for m in methods])
text+='<p>Each search uses the same 16-vCPU worker and existing reference cache. Allocation includes compute, provisioned EBS and public IPv4 during search wall time; it excludes staging, analysis, archival and idle time and is already included in lifetime spend. The reference occupies 318,551,974,004 bytes. No new persistent search index was built. Peak temporary disk use was not measured. These are sequential warm-reference observations.</p>'
text+='<h2>Per-mate natural evidence at E≤10⁻⁵</h2>'
groups={m:tax['natural_status_counts'][m+'/strict_E1e-5'] for m in methods}
columns=sorted({k for x in groups.values() for k in x})
text+=table(['Method']+columns,[[m]+[groups[m].get(k,0) for k in columns] for m in methods])
text+='<p>Natural evidence distinguishes artificial constructs from cellular or other nonviral matches. Viral-only describes returned eligible near-best taxa, not proven infection or a correct species assignment. Unknown taxonomy remains explicit.</p>'
changes={k:v for k,v in tax['natural_LCA_comparisons'].items() if k.endswith('/strict_E1e-5')}
columns=sorted({k for x in changes.values() for k in x})
text+=table(['Comparison']+columns,[[k]+[v.get(c,0) for c in columns] for k,v in sorted(changes.items())])
text+='<p>The unchanged count refers to LCA alone. Natural-status changes are counted separately and can overlap the LCA categories; a mate can change natural status while retaining the same LCA. Broader or narrower taxonomy is a description of resolution, not an accuracy label.</p>'
text+='<h2>Pair evidence by frozen stratum</h2>'
groups={k:v for k,v in pair['groups'].items() if '/strict_E1e-5/' in k or k.startswith('protein_NR/')}
text+=counts_table(groups,'status')
text+='<details><summary>Artificial-construct and natural-evidence separation</summary>'+counts_table(groups,'natural_status')+'</details>'
text+='<p>NR rows reuse the prior completed protein search and exactly reproduce all 232 frozen pair verdicts. They are contextual evidence: NR is a different reference and retains its original minimum 25-amino-acid filter and 90% near-best rule. The core_nt comparisons use best HSP per subject, at most 10 dense score ranks and 90% near-best. NR and core_nt costs or verdicts must not be interpreted as a controlled reference-equivalent comparison.</p>'
text+='<p>A separate reproducibility check on 746 diagnostic pairs changed 47 viral/nonviral statuses and 74 taxonomic assignments when the same reads were searched in a smaller batch with an automatically selected alternative DIAMOND mode. Those fractions cannot be extrapolated to this panel. The frozen pooled NR evidence is retained here with that limitation.</p>'
text+='<h2>Pairs recovered by translation</h2>'
gains=transitions['newly_viral_only_pairs']
text+=table(['Pair','Translated natural assignment','Rank','Aligned mates','NR status'],[[x['pair'],x['translated_name'],x['translated_rank'],x['translated_aligned_mates'],x['NR_status']] for x in gains])
losses=transitions['raw_viral_only_losses']
assert all(x['right_natural_status']=='viral_only' and x['translated_artificial_taxa']>0 and x['translated_cellular_taxa']==0 and x['translated_unknown_taxa']==0 for x in losses)
text+=f'<p>These are all pairs changing from no strict nucleotide alignment to natural viral-only evidence. The selected panel has one pair per sample, so these observations cannot establish new three-exemplar calls. All {len(losses)} raw viral-only-to-mixed transitions across both strata remain natural viral-only: their additional nonviral competitors are artificial constructs.</p>'
text+=table(['Mate','Aligned amino acids','AA identity %','Nucleotide identity %','Synonymous codon changes','Longest exact nucleotide stretch'],[[h['query'],h['alignment_aa'],h['reported_aa_identity_pct'],f"{h['aligned_nt_identity_pct']:.2f}",h['synonymous_changed_codons'],h['longest_exact_nt_run_in_alignment']] for x in gains for h in x['reconstructed_mates']])
text+='<p>Both synonymous changes and amino-acid substitutions occur in these alignments. The nucleotide-stretch measurement applies only within the reconstructed alignment; it does not identify the sole reason another method missed the read. '+link(d/'transition-analysis/summary.json','Complete transition and recovered-pair evidence')+'</p>'
text+='<h2>Translated sequence reconstruction</h2>'
rows=audit['rows']
text+=f'<p>For {len(rows)} changed mates with a strict translated hit, one deterministic highest-score subject was reconstructed from its original nucleotide coordinates. {sum(x["synonymous_changed_codons"]>0 for x in rows)} reconstructed alignments contain synonymous codon changes. This explains sequence-level differences; it does not prove the taxonomic verdict.</p>'
text+='<p>'+link(a.codon_audit,'Every reconstructed alignment, source DNA, frames, and extraction command')+'</p>'
text+='<h2>Selective translated checking</h2>'
selected=[x for x in routing['counts'] if x['threshold']=='strict_E1e-5' and x['comparison']=='natural_verdict']
text+=table(['Routing rule','Reads routed','Changed natural verdicts captured','Changes missed'],[[x['rule'],x['routed_queries'],x['captured_changes'],x['missed_changes']] for x in selected])
text+='<p>This retrospective simulation chooses between already measured per-read verdicts. It does not measure a faster search or establish which verdict is correct. Rules were chosen after the small pilot. '+link(d/'selective-validation/summary.json','Routing definitions, all counts and limitations')+'</p>'
text+='<h2>Reporting limits and exact provenance</h2>'
text+=table(['Method','Queries reaching 250 returned subjects','Maximum returned subjects'],[[m,tax['subject_cap_diagnostic'][m]['queries_with_at_least250_returned_subjects'],tax['subject_cap_diagnostic'][m]['maximum_returned_subjects']] for m in methods])
text+='<p>The subject cap can hide competing references. E≤10⁻⁵ is a postfilter of returned E≤10 hits, not a separately executed search. The input consists of separate original mates; production mate merging and complete end-to-end call generation were not replayed here. Different programs’ significance scores are not assumed interchangeable.</p>'
for path,label in [(d/'input-manifest.json','Frozen query selection and hashes'),(d/'analysis/per-mate.tsv','Every mate verdict at both thresholds'),(d/'analysis/disagreements.tsv','Every mate disagreement'),(d/'pair-analysis/per-pair.tsv','Every pair verdict'),(d/'pair-analysis/disagreements.tsv','Every pair disagreement'),(a.costs,'Completed-search commands, timings, source commits and cost accounting')]:
    text+='<p>'+link(path,label)+'</p>'
for m in methods:text+='<p>'+link(d/m/'protocol.json',m+' command and source commit')+'</p>'
text+='</html>'
with a.output.open('x') as f:f.write(text)
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':inputs,'report_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'completed_searches_required':methods}
with a.output.with_suffix('.provenance.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(a.output)

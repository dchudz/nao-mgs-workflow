"""Render the completed anchor analysis, keeping evidence and both policies explicit."""
import csv,html,json,subprocess
from pathlib import Path
r=Path(__file__).resolve().parent;d=r/'anchor-analysis/result-v1';s=json.loads((d/'summary.json').read_text());rows=list(csv.DictReader((d/'species-calls.tsv').open(),delimiter='\t'))
esc=lambda x:html.escape(str(x))
def table(headers,rows):return '<table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'
counts=s['counts'];raw=counts['nt/full/raw'];nat=counts['nt/full/exclude_constructs'];sh=counts['nt/first5m/raw']
statusrows=[]
for key in ['viral_only','mixed_viral_nonviral','nonviral_only','no_alignment','unknown_taxonomy']:
 statusrows.append([key,raw['statuses'].get(key,0),nat['statuses'].get(key,0),sh['statuses'].get(key,0)])
main=[x for x in rows if x['arm']=='nt' and x['depth']=='full' and x['policy']=='raw' and x['protein_call']=='True']
new=[x for x in main if x['gained_protein_call']=='True']
calltable=table(['Species (taxid)','Human annotation','Native clusters','Protein clusters','Protein-only clusters','Native comparison'],[[f"{x['name']} ({x['species']})",x['human_annotation'],x['native_clusters'],x['protein_clusters'],x['protein_only_clusters'],'new' if x['gained_protein_call']=='True' else 'already called'] for x in main])
depthrows=[]
for depth in ['first5m','full']:
 for policy in ['raw','exclude_constructs']:
  z=counts[f'nt/{depth}/{policy}'];c=[x for x in rows if x['arm']=='nt' and x['depth']==depth and x['policy']==policy and x['protein_call']=='True'];g=[x for x in c if x['gained_protein_call']=='True']
  depthrows.append([depth,policy,z['candidate_pairs'],z['protein_only_pairs'],z['protein_only_viral_only'],len(c),len(g),sum(x['human_annotation']=='1' for x in g)])
body=f'''<!doctype html><html><head><meta charset="utf-8"><title>Full-depth anchor NR evidence</title><style>body{{font:16px/1.5 system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#15202b}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:9px;border:1px solid #ccd5df;text-align:left}}th{{background:#eaf0f7}}h2{{margin-top:32px}}code{{word-break:break-all}}.note{{background:#f0f4f8;padding:16px}}</style></head><body>
<h1>Full-depth anchor: SRR39694979</h1><p>Analysis of preserved search results, 23 October 2026. No new sequence searches or cloud workers.</p>
<p><strong>{raw['protein_only_viral_only']:,} of the {raw['protein_only_pairs']:,} protein-only exact pairs have viral-only NR evidence.</strong> Excluding artificial constructs gives {nat['protein_only_viral_only']:,}. The strict species-concordant, both-NR-mates, three-similarity-cluster criterion supports {len(main)} protein species calls, of which {len(new)} are new relative to native in this sample.</p>
<h2>NR evidence for all 22,359 candidates</h2>{table(['Evidence class','Full depth: original','Full depth: exclude constructs','First 5 million: original'],statusrows)}
<p>Full depth contains 247,633,776 input pairs. The original 817,744-query NR search includes this anchor and logged <code>Algorithm: Double-indexed</code>. This analysis uses that frozen pooled output, not an auto-selected standalone rerun. Per-query near-best taxon sets reconstruct all {s['raw_NR_parity_pairs']:,} saved anchor classifications exactly.</p>
<p>Construct exclusion removes taxid 81077 and descendants from the already retained near-best taxa. It does not revisit lower-scoring alignments or reinterpret unknown taxa. Of all candidates, 5,435 have only viral and artificial-construct competition; every one was already present in native. Original classifications remain available alongside this separate policy.</p>
<h2>Species calls at full depth: original NR policy</h2>{calltable}
<p>Human annotation is retained from the frozen reference: 1 = match, 0 = inconsistent with human infection, 2 = unclear; NA = unavailable. None establishes infection in this sample. Protein-only clusters count supports containing eligible protein-only pairs; all protein clusters determine the protein-call threshold. Native supports are not symmetrically NR-validated.</p>
<h2>Depth comparison</h2>{table(['Depth','NR policy','Candidates','Protein-only pairs','Protein-only viral-only','Protein species calls','New vs native','Human-annotated new'],depthrows)}
<p>The shallow comparison selects original SRA spot IDs 1–5,000,000 from the full-depth evidence and reclusters those reads. It is a post-hoc within-sample depth comparison, not a separately executed shallow pipeline. The anchor was excluded from the 1,000-sample heldout list, whose four blocks were distributed across each run; this consecutive prefix is not that sampling design. It can show which retained evidence is absent at lower depth, not prove equivalent shallow-run behavior.</p>
<h2>Method and reproducibility</h2><p>Both NR mates must align and have combined viral-only near-best evidence. NR and viral-reference LCAs must resolve to the same species. The original pinned <code>nao_dedup</code> commit <code>404c7080e4f5a73d0fbdd02e535858c74b0ebeb0</code> clusters a joint pool of native retained pairs and all viral-only protein candidates, before mate/species restrictions, using real recovered qualities and deterministic sequence-hash order. Each depth and policy forms its own pool. The local implementation reproduced the previously saved worker fixture exactly.</p>
<p>This is the same supplemental similarity-only criterion as the 200-sample analysis. It does not reproduce the production alignment-coordinate duplicate stage, establish biological truth, or constitute an end-to-end pipeline comparison. Each read keeps its own validation verdict; no verdict is transferred through a cluster.</p>
<p>NT20 sensitivity-baseline counts and all below-threshold species are in the complete table. Runtime: {s['seconds']:.1f} local seconds. Source SHA-256: <code>{s['source_sha256']}</code>.</p>
<ul><li><a href="anchor-analysis/result-v1/summary.json">Summary, input hashes, commands, and limits</a></li><li><a href="anchor-analysis/result-v1/species-calls.tsv">All species, depths, policies, and native arms</a></li><li><a href="anchor-analysis/result-v1/pair-evidence.tsv">Every anchor candidate and classification</a></li><li><a href="anchor-analysis/independent-saved-table-counts.json">Independent check from prior saved tables</a></li></ul>
</body></html>'''
(r/'anchor-nr-result.html').write_text(body)
print(json.dumps({'full_raw_species_calls':len(main),'full_raw_new_species_calls':len(new),'new':new,'report':str(r/'anchor-nr-result.html')}))

"""Render completed same-reference measurements without inferring biological truth."""
import argparse
import csv
import datetime
import hashlib
import html
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
r = parser.parse_args().root.resolve()
p = r / 'pilot24'
analysis = p / 'analysis'
summary = json.loads((analysis / 'summary.json').read_text())
taxonomy = json.loads((analysis / 'taxonomy-description/summary.json').read_text())
audit = json.loads((analysis / 'codon-audit.json').read_text())
manifest = json.loads((p / 'input-manifest.json').read_text())
rate = 3.3345 + (2048 * .08 + 3000 * .005 + 875 * .04) / 720 + .005
methods = ['megablast', 'blastn', 'tblastx']
metrics, protocols = [], {}
for method in methods:
    result = summary['resources'][method]
    assert result['status'] == 'complete'
    protocols[method] = json.loads((p / method / 'protocol.json').read_text())
    counts = taxonomy['natural_status_counts'][method + '/strict_E1e-5']
    metrics.append(dict(method=method, query_mates=summary['mates'], wall_seconds=result['wall_seconds'], CPU_hours=(result['user_seconds'] + result['system_seconds']) / 3600, modeled_allocation_USD=result['wall_seconds'] / 3600 * rate, natural_viral_only_mates=counts.get('viral_only', 0), no_alignment_mates=counts.get('no_alignment', 0), mates_with_artificial_hits=taxonomy['rows_with_artificial_hits'][method + '/strict_E1e-5'], mates_at250_subject_cap=taxonomy['subject_cap_diagnostic'][method]['queries_with_at_least250_returned_subjects'], raw_alignment_rows=result['raw_alignment_rows'], search_source_commit=protocols[method]['source_commit']))
changes = [x for x in csv.DictReader((analysis / 'taxonomy-description/natural-LCA-changes.tsv').open(), delimiter='\t') if x['left'] == 'blastn' and x['right'] == 'tblastx' and x['threshold'] == 'strict_E1e-5']
seqs, q = {}, None
for line in (p / 'queries.fna').read_text().splitlines():
    if line.startswith('>'):
        q = line[1:].split()[0]
        seqs[q] = ''
    else:
        seqs[q] += line
lost = {x['query'] for x in changes if x['relation'] == 'lost_alignment'}
losses = []
for q in sorted(lost):
    best = {}
    for method in methods:
        with (p / method / 'alignments.tsv').open() as f:
            rows = [x for x in csv.reader(f, delimiter='\t') if x[0] == q]
        best[method] = max(rows, key=lambda x: float(x[3])) if rows else None
    losses.append(dict(query=q, length_nt=len(seqs[q]), maximum_full_frame_amino_acids=len(seqs[q]) // 3, best_returned_rows=best))
result = dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), pairs=summary['pairs'], mates=summary['mates'], metrics=metrics, strict_blastn_to_tblastx=taxonomy['natural_LCA_comparisons']['blastn -> tblastx/strict_E1e-5'], codon_audited_mates=len(audit['rows']), audited_mates_with_synonymous_codon_changes=sum(x['synonymous_changed_codons'] > 0 for x in audit['rows']), short_read_losses=losses, NR_cost_comparison='Pending completed NR search; not inferred from partial progress', raw_evidence_private_snapshot='s3://pipeline-lab-839993661825-us-east-1/runs/translated-core-nt-20261022/pilot-snapshots/original-three-raw/')
assert (result['pairs'], result['mates'], result['codon_audited_mates'], result['audited_mates_with_synonymous_codon_changes']) == (24, 43, 15, 8)
assert result['strict_blastn_to_tblastx'] == dict(unchanged=32, natural_status_changed=2, more_specific_LCA=4, broader_LCA=5, lost_alignment=2)
assert len(losses) == 2 and len({x['query'].rsplit('/', 1)[0] for x in losses}) == 1
assert seqs[losses[0]['query']].translate(str.maketrans('ACGT', 'TGCA'))[::-1] == seqs[losses[1]['query']]
for x in losses:
    assert x['length_nt'] == 68
    assert (int(x['best_returned_rows']['tblastx'][5]), float(x['best_returned_rows']['tblastx'][4])) == (21, .002)
    assert (int(x['best_returned_rows']['blastn'][5]), float(x['best_returned_rows']['blastn'][6]), float(x['best_returned_rows']['blastn'][4])) == (64, 100, 2.57e-22)
(r / 'pilot-result-summary.json').write_text(json.dumps(result, indent=2) + '\n')
e = lambda value: html.escape(str(value))
def table(headers, rows):
    return '<table><thead><tr>' + ''.join('<th>' + e(x) + '</th>' for x in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>' + e(x) + '</td>' for x in row) + '</tr>' for row in rows) + '</tbody></table>'
cost_table = table(['Method', 'Wall minutes', 'CPU hours', 'Modeled allocation', 'Natural hits: viral only', 'No strict hit', 'Artificial matches', 'At 250-subject cap'], [[x['method'], f"{x['wall_seconds'] / 60:.2f}", f"{x['CPU_hours']:.2f}", f"${x['modeled_allocation_USD']:.3f}", x['natural_viral_only_mates'], x['no_alignment_mates'], x['mates_with_artificial_hits'], x['mates_at250_subject_cap']] for x in metrics])
change_table = table(['Read mate', 'Stratum', 'Sensitive BLASTN', 'TBLASTX', 'Relationship'], [[x['query'], x['stratum'], x['left_name'], x['right_name'] or 'No strict hit', x['relation']] for x in changes])
history = table(['Method', 'Search source commit', 'Raw rows'], [[x['method'], x['search_source_commit'], x['raw_alignment_rows']] for x in metrics])
doc = f'''<!doctype html><html><head><meta charset="utf-8"><title>Same-reference translated pilot</title><style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1150px;margin:42px auto;padding:0 24px;color:#17212b}}h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0}}td,th{{padding:9px;border-bottom:1px solid #ccd4dd;text-align:left;overflow-wrap:anywhere}}th{{background:#eef3f8}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}}.note{{background:#f2f5f8;padding:15px}}a{{color:#155ba7}}</style></head><body>
<h1>Same-reference translated pilot</h1>
<p>Completed 43 mate queries from 24 pairs and 24 samples, selected before search outcomes: eight pairs each from baseline controls, amino-acid-screen additions, and NT20 additions. Every method searched the same complete core_nt snapshot. This compares validation of already selected reads; it does not measure a new detection gain.</p>
<h2>Measured searches</h2>{cost_table}
<p>All searches used 16 threads on the same worker. Allocation is wall time times${rate:.4f}/hour for EC2, configured EBS and IPv4, using a 720-hour month. These amounts are already included in the lifetime budget ledger. They exclude setup, failed earlier attempts, analyses and retained S3 storage. A shared existing reference cache was used; brief profiling and tiny characterization jobs overlapped TBLASTX.</p>
<p>Taxonomic columns use E≤10⁻⁵ as a postfilter of the original E≤10 output, the top 10 distinct score ranks, and 90% of the best score. They describe returned hits. The source pipeline’s LCA calculation is unchanged. Artificial constructs are counted separately.</p>
<h2>Natural assignment changes</h2>
<p>Relative to sensitive BLASTN, TBLASTX retained 32 natural LCAs, returned four more specific LCAs and five broader LCAs, and lost strict hits for two mates of one pair. Six mates additionally matched artificial constructs; those matches are not evidence of natural bacterial or human competitors.</p>
{change_table}
<h2>Sequence checks</h2>
<p>Fifteen changed mates with retained strict TBLASTX hits were checked against extracted original reference DNA. Eight had synonymous codon differences in the selected highest-scoring alignment. This verifies observed coding differences; it does not prove that any assignment is correct or that a synonymous change caused a nucleotide-search miss.</p>
<p>The two strict-hit losses are 68 nt mates, allowing at most 22 amino acids per full frame. Their best returned TBLASTX hits align 21 amino acids with E=0.002. Sensitive BLASTN instead finds a perfect 64 nt alignment with E=2.57×10⁻²². The two mate sequences are reverse complements, so these are not two independent examples.</p>
<h2>Limits of these measurements</h2>
<p class="note">The 250-subject reporting cap was reached for 38 of 43 TBLASTX queries and 27 of 43 sensitive BLASTN queries. Omitted competing references can affect LCA resolution. This pilot does not establish complete competitor enumeration, known-truth accuracy, infection, or population sensitivity. NR cost-per-query comparison remains pending its completed search.</p>
<p>Translated alignment lengths and identities use amino acids; nucleotide methods use bases. The characterized translated coverage filter uses peptide-frame length, and raw TBLASTX coverage fields can differ from nucleotide span. Both frame and nucleotide coverage are preserved in the alignment tables. Standard genetic code 1 was used.</p>
<h2>Evidence and history</h2>{history}
<p><a href="pilot24/analysis/per-mate.tsv">Original per-mate results</a> · <a href="pilot24/analysis/taxonomy-description/per-mate.tsv">Artificial/natural distinction</a> · <a href="pilot24/analysis/codon-audit.json">Source-DNA codon audit</a> · <a href="pilot-result-summary.json">Machine-readable summary</a></p>
<details><summary>Frozen inputs and exact commands</summary><pre>{e(json.dumps(manifest, indent=2))}</pre><pre>{e(json.dumps({k:v['command'] for k,v in protocols.items()}, indent=2))}</pre></details>
<p>Generated {e(result['generated_utc'])}. Private raw snapshot: {e(result['raw_evidence_private_snapshot'])}</p></body></html>'''
(r / 'same-reference-pilot.html').write_text(doc)
print(json.dumps({k: v for k, v in result.items() if k not in ['metrics', 'short_read_losses']}))

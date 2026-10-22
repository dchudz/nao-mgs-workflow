"""Render completed matched same-reference searches without conflating scoring and retrieval."""
import argparse,csv,hashlib,html,json,re,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('pilot',type=Path);p.add_argument('output',type=Path);p.add_argument('--analysis',default='analysis-query-fast');a=p.parse_args()
def load(path): return json.loads(path.read_text())
def wall(d):
 s=(d/'time.txt').read_text(); v=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)',s)[1];parts=list(map(float,v.split(':')));return sum(x*60**i for i,x in enumerate(parts[::-1]))
def esc(x):return html.escape(str(x))
def link(path,text):return '<a href="'+esc(path)+'">'+esc(text)+'</a>'
tax=load(a.pilot/a.analysis/'taxonomy-description/summary.json')
methods=['megablast','blastn','tblastx','query-tblastn-fast']
if (a.pilot/'tblastx-word4/complete.json').exists() and 'tblastx-word4/strict_E1e-5' in tax['natural_status_counts']:methods.append('tblastx-word4')
rows=[];commands=[];facts=[]
for method in methods:
 d=a.pilot/method; proto=load(d/'protocol.json'); done=load(d/'complete.json'); assert hashlib.sha256((d/'alignments.tsv').read_bytes()).hexdigest()==done['output_sha256']
 seconds=wall(d);count=tax['natural_status_counts'][method+'/strict_E1e-5'].get('viral_only',0);caps=tax['subject_cap_diagnostic'][method]['queries_with_at_least250_returned_subjects']
 rows.append(f'<tr><td>{esc(method)}</td><td>{seconds:,.2f}</td><td>${seconds/3600*3.6365:.3f}</td><td>{count}/43</td><td>{caps}/43</td></tr>')
 commands.append('<details><summary>'+esc(method)+' command and source commit</summary><p>'+esc(proto['source_commit'])+'</p><pre>'+esc(json.dumps(proto['command'],indent=2))+'</pre></details>')
 facts.append({'method':method,'wall_seconds':seconds,'allocated_compute_EBS_IPv4_USD':seconds/3600*3.6365,'natural_viral_mates_strict':count,'query_count':43,'subject_cap_queries':caps,'source_commit':proto['source_commit'],'output_sha256':done['output_sha256']})
shared=load(a.pilot/'shared-HSP-comparison/summary.json')
htmltext='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Faster translated validation on the same reference</title><style>body{font:17px/1.55 system-ui;max-width:1100px;margin:45px auto;padding:0 24px;color:#172332}table{border-collapse:collapse;width:100%;font-size:15px}td,th{border:1px solid #cad2dc;padding:9px;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}h1{line-height:1.15}a{color:#165d9c}.note{background:#f0f4f8;padding:18px}</style>
<h1>Faster translated validation on the same reference</h1>
<p>The six-frame query translation followed by <code>tblastn-fast</code> completed in 2,095.23 seconds, versus 11,408 seconds for default TBLASTX: a 5.44-fold reduction in wall time on the same 43 nucleotide reads and the same comprehensive <code>core_nt</code> reference. It also changes reported significance on identical alignments.</p>
<p class="note">This is a validation-method pilot on 24 selected pairs (43 mates), not a new detection run or an estimate of population sensitivity. The queries were frozen before these search outcomes. They include baseline, amino-acid-screen additions and NT20 additions. There are no known-truth labels.</p>
<h2>Matched measured searches</h2><table><tr><th>Method</th><th>Search seconds</th><th>Allocated cost</th><th>Natural viral mates, E ≤ 10⁻⁵</th><th>250-subject cap reached</th></tr>'''+''.join(rows)+'''</table>
<p>Each run used 16 threads on the same x2iedn.4xlarge worker, with the same pinned BLAST 2.17 container and unchanged full nucleotide database: 131,855,990 sequences and 1,088,816,606,629 nucleotides, 318,551,974,004 bytes in the frozen file inventory. No new reference index was built. Allocation uses $3.6365/hour for compute, provisioned EBS and IPv4; S3 and transfer are excluded. These allocations are already included in the lifetime budget ledger.</p>
<p>All methods searched with E ≤ 10 and were then filtered at E ≤ 10⁻⁵. Maximum 250 subjects and five HSPs per subject, with 30% query coverage. The fast method translates all six complete frames, including stop characters, and pools output back into original nucleotide-read identifiers and coordinates. Per-frame coverage and frame pooling are not identical to native TBLASTX behavior. Artificial constructs are separated from natural cellular competitors in the table; original unmodified classifications remain in the evidence.</p>
<h2>Retrieval and reported significance differ</h2>
<p>Of 28,285 byte-identical shared alignments (same read, reference accession, coordinates, frames and both aligned amino-acid strings), 22,786 have different bit scores. Scores in the fast method are between 0 and 32 bits higher. A total of 1,324 identical alignments cross the strict E-value cutoff; 1,180 contain a stop character. This comparison does not establish that one set of reported statistics is better calibrated.</p>
<p>The two apparent strict-cutoff gains over TBLASTX are mates of one 68-nucleotide pair, <code>SRR39695366.581592828</code>. Both methods already retrieve the same 21-amino-acid alignment to <code>PQ510853.1</code>: <code>VSRGY*SNTYNCD*EEH*FNT</code>. TBLASTX reports 51.9 bits and E=0.002; the fast method reports 61.7 bits and E=1.13×10⁻⁶. These two threshold crossings are not newly retrieved sequence evidence. Default TBLASTX also retrieves a stop-free 21-amino-acid alignment for this pair with E=0.002.</p>
<p>Among selected near-best alignments, 22 of 43 fast-method reads have only stop-containing alignments, compared with five of 41 strict-passing TBLASTX reads. This describes the retained near-best set, not an absence of stop-free alignments elsewhere in the raw output. No stop filter was added after observing this result.</p>
<h2>Taxonomic disagreements</h2>
<p>Compared with nucleotide <code>blastn -task blastn</code>, the fast method gives the same natural taxonomic ancestor for 38 of 43 reads, a more specific ancestor for two, and a broader ancestor for three. Compared with default TBLASTX, 37 ancestors are unchanged, two are more specific, two are broader, and the two short mates cross the significance cutoff. Greater specificity is not independently verified accuracy. Frequent output-cap saturation means omitted competing references remain a limitation.</p>
<p>The original full TBLASTX pilot is not being scaled. A bounded same-query native TBLASTX pilot using word size four and threshold 18 tests whether a longer seed can reduce work while preserving scoring on identical retrieved alignments.</p>
<h2>Evidence and reproducibility</h2>'''
for path,label in [(a.pilot/a.analysis/'taxonomy-description/per-mate.tsv','Per-read natural classifications'),(a.pilot/a.analysis/'taxonomy-description/natural-LCA-changes.tsv','Every changed natural taxonomic assignment'),(a.pilot/'shared-HSP-comparison/shared-identical-alignments.tsv.gz','Identical alignments and their scores'),(a.pilot/'shared-HSP-comparison/per-query.tsv','Per-read retrieval and score comparison'),(a.pilot/a.analysis/'codon-audit-query-tblastn-fast-vs-tblastx.json','Changed-read source DNA and codon reconstruction')]:
 assert path.exists(),path
 htmltext+='<p>'+link(os.path.relpath(path.resolve(),a.output.parent.resolve()),label)+'</p>'
htmltext+=''.join(commands)+'</html>'
a.output.write_text(htmltext)
a.output.with_suffix('.json').write_text(json.dumps({'methods':facts,'limits':['Selected small panel, not population sensitivity.','Protein and nucleotide E-values are not interchangeable.','Costs are allocations, not an invoice.','Per-frame coverage and output pooling differ.','No new reference index was built.']},indent=2)+'\n')
print(a.output)

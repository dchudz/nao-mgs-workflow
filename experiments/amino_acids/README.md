# Amino-acid detection experiment

These are the private experimental workflows and analysis sources used for the
2026-10-21 comparison. They are not enabled by INDEX, RUN, or DOWNSTREAM.
The frozen protocol defines the primary comparisons and evidence thresholds.
Resource and secondary-analysis changes are recorded separately with the results.

The first comparison screens the same raw pairs with NT24, NT22, NT20, AA12,
and NT24 OR AA12. `postscreen.nf` retains the original nucleotide processing after
the initial screen. The second comparison uses the AA12 candidates, original
FASTP and contaminant filters, exact complete-pair deduplication, translated
protein alignment, taxonomic LCA, and broad protein-reference competition.
Protein evidence has its own output schema; it is not converted into nucleotide
alignment scores or coordinates.

`protein_queries.py` retains every original ID in an exact-pair membership table.
`protein-viral-summary.py` summarizes the fixed viral-reference alignments.
`protein-nr-summary.py` reports broad-reference viral, nonviral, mixed, unknown,
and unaligned evidence. `compare-protein.py` and `analyze-protein200.py` compare
complete original trimmed pairs, avoiding missing-mate fields in native TSVs.

The workflows reuse the repository's original modules. The single-worker
configuration uses local references, pinned containers, and disables Wave,
Fusion, and Tower. Task containers have no network; the host stages public
inputs and private results. Python analysis scripts retain the deployed worker
root `/home/ec2-user/experiment` for exact reproducibility. They are experimental
entry points, not a portable public interface.

The local experiment archive contains the deployment scripts, input manifests,
reference and executable hashes, task traces, all protocol addenda, results,
and resource records. Its path is
`/Users/david/Documents/ChatGPT/pipeline/.tmp/aa-experiment-20261021`.
The AWS copy is private. No results or code have been published.

Results remain in progress. Environmental assignments do not provide known-truth
sensitivity or specificity; added reads or species assignments do not establish
new infections or novel viruses. The final report must state actual completed
sample counts, read depth, broad-reference validation coverage, and resource use.

## Secondary repeat-filter experiment

An observed candidate-volume outlier motivated omitting AA12 words containing
ten or more copies of one residue. This removes 6,900 of 93,380,138 index words.
The initial diagnostic uses the five development samples and the triggering
outlier. A separate protocol fixes evaluation on the last 500 original manifest
samples before any of those samples started. Neither changes the primary
1,000-sample comparison or 200-sample protein path. Report this chronology and
all assignment losses; preservation is measured, not assumed.

The resource addenda record byte-identical full-depth aggregation with cyclic
garbage collection disabled and equivalence of bounded NR taxon-score storage.
The original NR summarizer is retained as an independent parity oracle.

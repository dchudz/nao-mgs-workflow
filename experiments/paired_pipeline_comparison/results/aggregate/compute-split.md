# Measured stage time

These are completed-unit wall times, not CPU times or a bill. Cold reference work, validation order, and sampled records can affect them. The 5M and 25M strata are nested and are kept separate.

| Pairs | Stage | Units | Median minutes | Range minutes |
|---|---|---:|---:|---:|
| 5M | aa-screen | 18 | 1.19 | 1.10–3.66 |
| 5M | baseline-run | 18 | 4.80 | 3.59–15.62 |
| 5M | prepare | 18 | 0.06 | 0.06–0.07 |
| 5M | production-downstream | 18 | 2.55 | 1.20–39.64 |
| 5M | protein-run | 18 | 3.30 | 2.69–12.72 |
| 5M | protein_variant-downstream | 18 | 5.93 | 1.82–19.11 |
| 5M | union | 18 | 0.00 | 0.00–0.01 |
| 25M | aa-screen | 9 | 3.89 | 3.71–4.01 |
| 25M | baseline-run | 9 | 4.25 | 2.58–5.85 |
| 25M | prepare | 9 | 0.09 | 0.06–0.15 |
| 25M | production-downstream | 9 | 3.70 | 1.42–7.65 |
| 25M | protein-run | 9 | 4.68 | 3.26–11.86 |
| 25M | protein_variant-downstream | 9 | 4.02 | 1.97–20.22 |
| 25M | union | 9 | 0.02 | 0.00–0.05 |

Acquisition of the entire public source archive, input sampling, boot, preservation, retries and idle periods are outside this table. Acquisition is shared by both arms and by nested depths. The execution-overhead receipts account for wrapper time outside the listed stages; they do not attribute that residual to a search.

The baseline run also produces QC/profiling outputs reused by the variant. Adding the protein stages to baseline time is the incremental cost of this paired experiment, not a standalone deployment cost ratio. A standalone variant would still need common QC/profiling work, whose cost is not isolated by these top-level stages.

Use runtime-context.tsv for measured overlaps with reference-storage checks and the disjoint worker-interval ledger for compute-cost bounds. No overlap duration is subtracted here.

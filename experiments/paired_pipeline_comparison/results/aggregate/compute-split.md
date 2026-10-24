# Measured stage time

These are completed-unit wall times, not CPU times or a bill. Cold reference work, validation order, and sampled records can affect them. The 5M and 25M strata are nested and are kept separate.

| Pairs | Stage | Units | Median minutes | Range minutes |
|---|---|---:|---:|---:|
| 5M | aa-screen | 14 | 1.19 | 1.10–3.66 |
| 5M | baseline-run | 14 | 4.80 | 3.59–15.62 |
| 5M | prepare | 14 | 0.06 | 0.06–0.07 |
| 5M | production-downstream | 14 | 2.55 | 1.20–39.64 |
| 5M | protein-run | 14 | 3.26 | 2.69–12.72 |
| 5M | protein_variant-downstream | 14 | 7.13 | 1.82–19.11 |
| 5M | union | 14 | 0.00 | 0.00–0.01 |
| 25M | aa-screen | 6 | 3.90 | 3.71–4.01 |
| 25M | baseline-run | 6 | 4.27 | 2.58–5.85 |
| 25M | prepare | 6 | 0.09 | 0.06–0.15 |
| 25M | production-downstream | 6 | 3.68 | 1.77–7.65 |
| 25M | protein-run | 6 | 5.66 | 3.48–11.86 |
| 25M | protein_variant-downstream | 6 | 5.64 | 1.97–20.22 |
| 25M | union | 6 | 0.01 | 0.00–0.05 |

Acquisition of the entire public source archive, input sampling, boot, preservation, retries and idle periods are outside this table. Acquisition is shared by both arms and by nested depths. The execution-overhead receipts account for wrapper time outside the listed stages; they do not attribute that residual to a search.

The baseline run also produces QC/profiling outputs reused by the variant. Adding the protein stages to baseline time is the incremental cost of this paired experiment, not a standalone deployment cost ratio. A standalone variant would still need common QC/profiling work, whose cost is not isolated by these top-level stages.

Use runtime-context.tsv for measured overlaps with reference-storage checks and the disjoint worker-interval ledger for compute-cost bounds. No overlap duration is subtracted here.

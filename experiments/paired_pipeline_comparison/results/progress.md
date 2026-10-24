# Frozen paired pipeline comparison

Report snapshot: 2026-10-24T08:29:34.713569+00:00.

8 of 12 feasibility comparisons are collected and accepted. The full frozen design contains 60 comparisons: 50 at 5 million pairs and 10 nested comparisons at 25 million pairs. A missing local report can still be running remotely.

Both arms use the same selected input pairs, source reference collection, production duplicate handling, validation mode and validation cap. The variant changes the upstream screen and alignment representation. This does not isolate alignment as the only cause of differences.

## Completed output inventory

These are final output records, not confirmed detections or species calls. “Only” compares pair identifiers between arms at the same depth. It does not compare shallow and deep identifiers.

| Run | Input pairs (million) | Production records | Protein records | Production only | Protein only | Stage minutes | Driver minutes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SRR38294967 | 5 | 2,861 | 3,408 | 8 | 555 | 37.9 | unrecorded |
| SRR39695348 | 5 | 112 | 158 | 0 | 46 | 74.5 | 74.6 |
| SRR39695348 | 25 | 594 | 834 | 1 | 241 | 16.3 | 27.9 |
| SRR38294785 | 5 | 344 | 620 | 0 | 276 | 22.2 | 22.3 |
| SRR37059981 | 5 | 931 | 1,396 | 0 | 465 | 16.0 | 16.1 |
| SRR38294137 | 5 | 691 | 1,639 | 1 | 949 | 20.3 | 20.4 |
| SRR39695239 | 5 | 2,986 | 4,182 | 6 | 1,202 | 65.0 | 65.3 |
| SRR39695239 | 25 | 14,746 | 20,877 | 29 | 6,160 | 23.4 | 25.3 |

Stage minutes omit preflight, inter-stage backups and driver overhead. Driver minutes include those, but omit acquisition and the final verified archive. Worker uptime captures all of them. Early replica timings overlap storage initialization; see [runtime context](aggregate/runtime-context.tsv). The same input QC/profiling work is shared between arms, so these are not standalone deployment timings.

## Interpretation

- The known-origin artificial control passed 24 positive pairs, one negative pair, and 12 expected duplicate groups using the production duplicate container.
- Every accepted unit has both arms complete, equal downstream parameters, output checks, and a checksum-verified private S3 archive.
- Bracken abundance is unavailable because the retained taxonomy/report combination has no requested D-rank rows. These are accepted assignment comparisons, not complete abundance reports.
- Clade counts precede validation. Validation annotates sampled records and does not turn all output records into confirmed positives.
- Standard validation agreement allows a parent label to agree with a more specific child label. It is not exact species agreement or accuracy.
- Output additions, losses and changes have no complete environmental truth labels. More output is not established to be better.
- Nested depth comparisons reuse observations from the same run. Differences of exemplar totals are not counts of new independent origins.

See [the interpretation record](../comparison-interpretation.md), [nested depth table](aggregate/nested-depth.tsv), and [full acceptance inventory](aggregate/runs.tsv). Individual run folders retain the standard comparison reports, source commit, commands, stage timings, acquisition receipts and private archive checksum.

## Remaining execution

4 feasibility comparisons remain uncollected. Expansion beyond this initial cohort is conditional on completed outputs and measured costs. No sample will be chosen or discarded because of its biological findings.

The existing budget is $500. The working allocation reserves $200 for new compute; posted historical charges, billing lag and noncompute retention are accounted for separately in [budget.json](../budget.json). Observed resource-time ledgers live in `../runtime-provenance/`. They are estimates and ceilings, not invoice reconciliation.

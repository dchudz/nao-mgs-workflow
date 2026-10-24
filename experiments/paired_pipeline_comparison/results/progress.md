# Frozen paired pipeline comparison

Report snapshot: 2026-10-24T12:22:36.351052+00:00.

15 of 60 planned comparisons are collected and accepted, including 12 of 12 feasibility comparisons. The full frozen design contains 60 comparisons: 50 at 5 million pairs and 10 nested comparisons at 25 million pairs. A missing local report can still be running remotely.

Both arms use the same selected input pairs, source reference collection, production duplicate handling, validation mode and validation cap. The variant changes the upstream screen and alignment representation. This does not isolate alignment as the only cause of differences.

## Completed output inventory

These are final output records, not confirmed detections or species calls. “Only” compares pair identifiers between arms at the same depth. It does not compare shallow and deep identifiers. Reassigned shared records occur in both arms with different labels; these are separate from added or lost records.

| Run | Input pairs (million) | Production records | Protein records | Production only | Protein only | Shared records reassigned | Stage minutes | Driver minutes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SRR38294967 | 5 | 2,861 | 3,408 | 8 | 555 | 1,652 | 37.9 | unrecorded |
| SRR39695144 | 5 | 9,871 | 11,048 | 87 | 1,264 | 6,660 | 21.8 | 22.0 |
| SRR39695348 | 5 | 112 | 158 | 0 | 46 | 84 | 74.5 | 74.6 |
| SRR39695348 | 25 | 594 | 834 | 1 | 241 | 438 | 16.3 | 27.9 |
| SRR38294785 | 5 | 344 | 620 | 0 | 276 | 216 | 22.2 | 22.3 |
| SRR37059981 | 5 | 931 | 1,396 | 0 | 465 | 556 | 16.0 | 16.1 |
| SRR38294137 | 5 | 691 | 1,639 | 1 | 949 | 418 | 20.3 | 20.4 |
| SRR38294805 | 5 | 2,974 | 3,893 | 2 | 921 | 2,182 | 19.7 | 20.3 |
| SRR38294805 | 25 | 15,063 | 19,542 | 21 | 4,500 | 11,169 | 27.5 | 29.4 |
| SRR37367647 | 5 | 26 | 43 | 0 | 17 | 20 | 15.5 | 15.7 |
| SRR39695239 | 5 | 2,986 | 4,182 | 6 | 1,202 | 2,111 | 65.0 | 65.3 |
| SRR39695239 | 25 | 14,746 | 20,877 | 29 | 6,160 | 10,327 | 23.4 | 25.3 |
| SRR38294858 | 5 | 1,482 | 2,008 | 9 | 535 | 1,013 | 21.4 | 21.8 |
| SRR38294103 | 5 | 1,222 | 1,761 | 1 | 540 | 833 | 18.5 | 18.7 |
| SRR38294653 | 5 | 29 | 388 | 0 | 359 | 14 | 29.1 | 29.2 |

Stage minutes omit preflight, inter-stage backups and driver overhead. Driver minutes include those, but omit acquisition and the final verified archive. Worker uptime captures all of them. Early replica timings overlap storage initialization; see [runtime context](aggregate/runtime-context.tsv). The same input QC/profiling work is shared between arms, so these are not standalone deployment timings.

## Interpretation

- The known-origin artificial control passed 24 positive pairs, one negative pair, and 12 expected duplicate groups using the production duplicate container.
- Every accepted unit has both arms complete, equal downstream parameters, output checks, and a checksum-verified private S3 archive.
- Bracken abundance is unavailable because the retained taxonomy/report combination has no requested D-rank rows. These are accepted assignment comparisons, not complete abundance reports.
- Clade counts precede validation. Validation annotates sampled records and does not turn all output records into confirmed positives.
- Standard validation agreement allows a parent label to agree with a more specific child label. It is not exact species agreement or accuracy.
- Output additions, losses and changes have no complete environmental truth labels. More output is not established to be better.
- Nested depth comparisons reuse observations from the same run. Differences of exemplar totals are not counts of new independent origins.

See [the interpretation record](../comparison-interpretation.md), [nested depth table](aggregate/nested-depth.tsv), [measured stage split](aggregate/compute-split.md), and [full acceptance inventory](aggregate/runs.tsv). Individual run folders retain the standard comparison reports, source commit, commands, stage timings, acquisition receipts and private archive checksum.

## Remaining execution

45 of the full cohort comparisons remain uncollected. The feasibility gate passed and the fixed expansion was accepted within the budget. No sample is chosen or discarded because of its biological findings. Uncollected does not mean a zero result.

The existing budget is $500. The working allocation reserves $200 for new compute; posted historical charges, billing lag and noncompute retention are accounted for separately in [budget.json](../budget.json). Observed resource-time ledgers live in `../runtime-provenance/`. They are estimates and ceilings, not invoice reconciliation.

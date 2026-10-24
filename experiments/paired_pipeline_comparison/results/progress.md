# Frozen paired pipeline comparison

Report snapshot: 2026-10-24T22:02:10.875522+00:00.

38 of 60 planned comparisons are collected and accepted, including 12 of 12 feasibility comparisons. The full frozen design contains 60 comparisons: 50 at 5 million pairs and 10 nested comparisons at 25 million pairs. A missing local report can still be running remotely.

Both arms use the same selected input pairs, source reference collection, production duplicate handling, validation mode and validation cap. The variant changes the upstream screen and alignment representation. This does not isolate alignment as the only cause of differences.

## Completed output inventory

These are final output records, not confirmed detections or species calls. “Only” compares pair identifiers between arms at the same depth. It does not compare shallow and deep identifiers. Reassigned shared records occur in both arms with different labels; these are separate from added or lost records.

| Run | Input pairs (million) | Production records | Protein records | Production only | Protein only | Shared records reassigned | Stage minutes | Driver minutes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SRR38294967 | 5 | 2,861 | 3,408 | 8 | 555 | 1,652 | 37.9 | unrecorded |
| SRR39695225 | 5 | 3,624 | 4,302 | 7 | 685 | 2,429 | 22.8 | 23.0 |
| SRR39695144 | 5 | 9,871 | 11,048 | 87 | 1,264 | 6,660 | 21.8 | 22.0 |
| SRR39695144 | 25 | 49,581 | 55,537 | 404 | 6,360 | 33,570 | 33.7 | 34.6 |
| SRR35971518 | 5 | 38 | 311 | 1 | 274 | 18 | 13.6 | 13.8 |
| SRR35971518 | 25 | 204 | 1,566 | 2 | 1,364 | 96 | 16.0 | 16.3 |
| SRR33038009 | 5 | 156 | 460 | 1 | 305 | 77 | 23.1 | 23.2 |
| SRR39695348 | 5 | 112 | 158 | 0 | 46 | 84 | 74.5 | 74.6 |
| SRR39695348 | 25 | 594 | 834 | 1 | 241 | 438 | 16.3 | 27.9 |
| SRR38294785 | 5 | 344 | 620 | 0 | 276 | 216 | 22.2 | 22.3 |
| SRR37059981 | 5 | 931 | 1,396 | 0 | 465 | 556 | 16.0 | 16.1 |
| SRR36876827 | 5 | 232 | 444 | 0 | 212 | 116 | 20.0 | 20.2 |
| SRR36876827 | 25 | 1,121 | 2,115 | 2 | 996 | 554 | 35.4 | 35.7 |
| SRR38294137 | 5 | 691 | 1,639 | 1 | 949 | 418 | 20.3 | 20.4 |
| SRR38294805 | 5 | 2,974 | 3,893 | 2 | 921 | 2,182 | 19.7 | 20.3 |
| SRR38294805 | 25 | 15,063 | 19,542 | 21 | 4,500 | 11,169 | 27.5 | 29.4 |
| SRR35946910 | 5 | 8,667 | 13,084 | 117 | 4,534 | 2,159 | 22.1 | 22.3 |
| SRR39695024 | 5 | 654 | 1,193 | 2 | 541 | 498 | 19.6 | 19.7 |
| SRR37367643 | 5 | 946 | 1,459 | 5 | 518 | 645 | 24.1 | 24.3 |
| SRR35987549 | 5 | 867 | 1,689 | 1 | 823 | 617 | 20.8 | 21.0 |
| SRR38294666 | 5 | 2,097 | 2,899 | 6 | 808 | 1,516 | 19.4 | 19.5 |
| SRR37367647 | 5 | 26 | 43 | 0 | 17 | 20 | 15.5 | 15.7 |
| SRR37092454 | 5 | 404 | 1,005 | 1 | 602 | 228 | 23.9 | 24.1 |
| SRR38294124 | 5 | 290 | 2,318 | 1 | 2,029 | 132 | 30.5 | 30.6 |
| SRR38294124 | 25 | 1,509 | 11,549 | 7 | 10,047 | 665 | 19.8 | 20.0 |
| SRR38294667 | 5 | 45 | 45 | 3 | 3 | 19 | 17.9 | 18.2 |
| SRR35987661 | 5 | 54 | 328 | 0 | 274 | 28 | 18.9 | 19.0 |
| SRR37367662 | 5 | 819 | 1,551 | 3 | 735 | 535 | 17.7 | 17.9 |
| SRR37367662 | 25 | 4,292 | 7,696 | 8 | 3,412 | 2,721 | 22.8 | 23.3 |
| SRR38294842 | 5 | 415 | 955 | 1 | 541 | 253 | 17.4 | 17.5 |
| SRR38294842 | 25 | 2,021 | 4,842 | 9 | 2,830 | 1,239 | 22.5 | 22.8 |
| SRR39695239 | 5 | 2,986 | 4,182 | 6 | 1,202 | 2,111 | 65.0 | 65.3 |
| SRR39695239 | 25 | 14,746 | 20,877 | 29 | 6,160 | 10,327 | 23.4 | 25.3 |
| SRR39695102 | 5 | 527 | 1,065 | 2 | 540 | 375 | 15.8 | 16.0 |
| SRR39695102 | 25 | 2,690 | 5,407 | 8 | 2,725 | 1,873 | 19.7 | 20.1 |
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
- Equal downstream settings do not imply matched upstream scoring thresholds or false-positive operating points.
- Nested depth comparisons reuse observations from the same run. Differences of exemplar totals are not counts of new independent origins.

See [the interpretation record](../comparison-interpretation.md), [nested depth table](aggregate/nested-depth.tsv), [measured stage split](aggregate/compute-split.md), and [full acceptance inventory](aggregate/runs.tsv). Individual run folders retain the standard comparison reports, source commit, commands, stage timings, acquisition receipts and private archive checksum.

## Remaining execution

22 of the full cohort comparisons remain uncollected. The feasibility gate passed and the fixed expansion was accepted within the budget. No sample is chosen or discarded because of its biological findings. Uncollected does not mean a zero result.

The existing budget is $500. The working allocation reserves $200 for new compute; posted historical charges, billing lag and noncompute retention are accounted for separately in [budget.json](../budget.json). Observed resource-time ledgers live in `../runtime-provenance/`. They are estimates and ceilings, not invoice reconciliation.

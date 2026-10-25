# Final paired comparison report

The report covers 50 public sequencing runs at 5 million read pairs and ten nested comparisons at 25 million pairs. All 60 comparisons completed. These are aggregate measurements of the frozen comparison, not sensitivity or specificity estimates.

`WRITEUP.md` is the reviewed snapshot of the shared Claude/Codex draft. The tab-separated tables beside it support its numerical claims. Code remains in the parent experiment directory; each provenance file records method and input hashes. Absolute input paths identify the original working environment.

`final-aggregate-report-inputs.verified.json` identifies the private, checksum-verified archive containing per-comparison count tables, validation-category aggregates and recovered runtime traces. Extracting it restores the relative `species-audit`, `validation-categories` and `trace-recovery` directories used by the reporting scripts. The separately preserved `species-transitions.tsv` output is identified by `species-transitions.verified.json`; it is omitted from the presentation snapshot, not from preservation.

The scientific source, input receipts and verified per-comparison archive locations are recorded under `../results/`. Completion and shutdown records are under `../runtime-provenance/completion/`. `cost-snapshot.json` closes the compute intervals; its disk accrual is an observation rather than a final storage bill. All six volumes remain retained.

The branch and this report remain local and unpushed. The secondary NR audit was not performed.

# Final paired comparison report

The report covers 50 public sequencing runs at 5 million read pairs and ten nested comparisons at 25 million pairs. All 60 comparisons completed. These are aggregate measurements of the frozen comparison, not sensitivity or specificity estimates.

`WRITEUP.md` is the reviewed snapshot of the shared Claude/Codex draft. The tab-separated tables beside it support its numerical claims. Code remains in the parent experiment directory; each provenance file records method and input hashes. Absolute input paths identify the original working environment.

`final-aggregate-report-inputs.verified.json` identifies the private, checksum-verified archive containing per-comparison count tables, validation-category aggregates and recovered runtime traces. Extracting it restores the relative `species-audit`, `validation-categories` and `trace-recovery` directories used by the reporting scripts. The separately preserved `species-transitions.tsv` output is identified by `species-transitions.verified.json`; it is omitted from the presentation snapshot, not from preservation.

The scientific source, input receipts and verified per-comparison archive locations are recorded under `../results/`. Completion and shutdown records are under `../runtime-provenance/completion/`. `cost-snapshot.json` closes the compute intervals; its disk accrual is an observation rather than a final storage bill. All six volumes remain retained.

The branch and this report remain local and unpushed. The secondary NR audit was not performed.

## Addendum, 2026-10-25

Two statements above were accurate when written and are now superseded.

**"The branch and this report remain local and unpushed."** The branch was pushed to the public
fork `dchudz/nao-mgs-workflow` on 2026-10-25 and opened as a draft pull request, as an archival
artifact of the research rather than a proposal to merge. The commits were scanned for
credentials before the push: no access keys, secret or session tokens, GitHub tokens, private
keys, presigned URL signatures or local paths appear in any of the 8,433 blobs or 442 commit
messages. The AWS account ID appears in 346 files and was published deliberately. All sample
data is public SRA.

**"All six volumes remain retained."** All EBS volumes have since been deleted, after each was
snapshotted and the snapshot confirmed `completed`. Seven snapshots now hold the preserved
state, and the one remaining worker instance was terminated. `cost-snapshot.json` and the
closure records describe the state at closure, not the state now;
`runtime-provenance/completion/post-cleanup-inventory.json` records the state after cleanup.

A correction to the closure cost model is recorded with those inventories: `FullSnapshotSizeInBytes`
is not a billed size, and three of the snapshots were restored from a common parent and share
blocks, so snapshot storage cannot be derived by summing full sizes.

# Storage after experiment closure

Codex independently observed the following after Claude reported cleanup under David's authorization. Codex performed read-only AWS inventory calls and did not perform the cleanup.

At 2026-10-25T06:23:20.551150+00:00, five formerly detached volumes were absent, the original volume remained attached to the stopped primary instance, and six named snapshots were completed. See `post-cleanup-inventory.json`.

At 2026-10-25T11:39:16.302959+00:00, there were zero EBS volumes and zero pending, running, stopping, or stopped instances in us-east-1 and us-west-2. All seven named snapshots were completed in the standard tier. See `final-cleanup-inventory.json` for their IDs, source volumes, timestamps and sizes. These observations supersede the retained-resource state at closure; they do not change the historical compute or study-window EBS estimates.

The fresh original-volume snapshot is `snap-001bd9d85d254dbdb`, started at 2026-10-25T10:49:14.079000+00:00 and completed at 2026-10-25T11:36:31.606000+00:00. The earlier `snap-0c13ccf01c762eb8a` captured the original volume before the paired experiment. A difference between their full snapshot sizes does not measure unique file contents, changed blocks, or data absent from S3. No file-level comparison was performed here.

Snapshot storage continues accruing. `FullSnapshotSizeInBytes` describes full recoverable data size, not incremental billed storage. Snapshot cost cannot be established by summing these values. See the [AWS Snapshot API](https://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_Snapshot.html) and [snapshot storage documentation](https://docs.aws.amazon.com/ebs/latest/userguide/how_snapshots_work.html). The restored environment volumes have shared snapshot ancestry; the two scratch volumes are not established as members of that lineage.

The earlier 168.27 GB S3 inventory covered only `runs/paired-comparison-20261024/`, not the full bucket. This follow-up did not independently verify Claude's later S3 deletions, inventories, or rate estimates. It also did not reconcile account-wide billing against the study-window model; those scopes and charge components differ.

No experiment was restarted, no sequence analysis was added, and nothing was pushed or published by Codex during this follow-up.

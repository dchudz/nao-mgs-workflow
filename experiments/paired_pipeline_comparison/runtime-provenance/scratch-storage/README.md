# Input storage recovery

The two replica root volumes were still initializing from their snapshots. AWS
reported normal I/O health and initialization progress of 31% and 34% at the
recorded observation. This is context, not proof that initialization explains
all observed latency.

A disposable 512 MiB random-byte file on spot1 took 7.68 seconds to write, flush
and hash, followed by 91.49 seconds to read and hash after a file-local cache
hint. The same procedure on a new empty scratch volume took 2.06 and 1.93
seconds. Checksums matched in both checks. Other I/O was active and the random
payloads differed, so this is an operational diagnostic, not a controlled
storage benchmark or a pipeline speedup claim.

Two encrypted 512 GiB gp3 scratch volumes now hold temporary input preparation.
Each has 3,000 IOPS and 500 MiB/s throughput. Their combined modeled cost is
$3.73 per 30-day-prorated day, or $4.23 in the higher throughput-rate scenario.
This uses the existing noncompute reserve. References, indexes and pipeline
outputs remain on their existing disks. `scratch-storage.json` in each worker
folder records the volume ID, filesystem UUID and persistent mount.

The unfinished SRR38294858 preparation on spot1 was interrupted after its
validator continued reading slowly. SIGTERM was handled without exiting, so the
same verified process was stopped with SIGKILL. No comparison stage had begun.
The original source archive remains in `input-attempts/SRR38294858-restored-disk`;
its logs and the interruption are preserved in the checksum-verified private
archive named in `spot1/storage-retry-evidence-backup.json`. The separate original
sequential checksum was interrupted and has no completed digest.

Before interruption, a new copy of the same public S3 object on fresh scratch
storage completed downloading in 147.74 seconds, hashing in 126.59 seconds and
the unchanged VDB validation in 106.86 seconds. The source object's recorded
last modification predates both attempts. Its checksum is recorded in
`spot1/fresh-storage-preflight.json`; the final acquisition receipt must match.
The retry uses a new log/shard name containing only the two unfinished runs,
so completed comparison archives are not overwritten. Pipeline source remains
`c77cdde08b03834a9c1f54f3d80b0c9668bf9a83`, with the same sample selection,
paired-read sampler, reference, search settings and validation policy.

Spot2's current SRR37367647 preparation recovered, passed validation and finished
sampling on its original disk; it was not interrupted. The independent checksum
in its folder must match the completed acquisition receipt. Its next input,
SRR38294103, and spot1's next input, SRR38294653, use persistent directory
symlinks into their fresh scratch volumes. The logical paths and frozen commands
are unchanged. Both source-storage interventions and the retry consume worker
time and remain part of experiment spending; none is subtracted as free time.

Temporary scratch and replica disks must remain available until their unique
results and selected inputs have verified private backups. The original primary
reference disk and snapshot are separate retained assets.

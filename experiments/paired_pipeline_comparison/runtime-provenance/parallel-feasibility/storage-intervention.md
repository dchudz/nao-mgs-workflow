# Restored-volume startup

The first replica profiling tasks were waiting on memory-mapped index reads.
The active mapping was verified as `/scratch/kraken_db_9f172f33/hash.k2d`,
backed by the retained experiment scratch directory. Low CPU utilization did
not mean idle workers.

A read-only checksum pass touched this same 117,958,174,752-byte file on all
three workers. The retained worker finished in 172.13 seconds. The restored
replicas initially progressed much more slowly. Bounded parallel reads were
then added on the two replicas, using 32 readers, 1 MiB blocks and disjoint
1 GiB regions. Files were opened read-only; size and modification time were
checked afterward. No reference, search setting, taxonomy or calling rule
changed. The two initialization passes took 129.82 and 490.27 seconds.

All three whole-file SHA-256 results match. The replica checksum elapsed times
include the initial serial read period and concurrent parallel initialization;
they are not a controlled serial-versus-parallel benchmark. Their first
pipeline executions overlap this storage intervention and must not be called
steady-state deployment timings. All worker time remains in the spending ledger.

AWS documents increased latency while restoring snapshot blocks and supports
read-only initialization: [Initialize Amazon EBS volumes](https://docs.aws.amazon.com/ebs/latest/userguide/initalize-volume.html).
The experiment reads only the needed file, not the whole device. The source
scripts, timestamps and checksums are retained alongside these records.

The retained 4,313,022,605-byte alignment index also matches its previously
recorded checksum on all workers. Reading it did not rebuild or alter it.

A separate full-cache check covered 836 validation files totaling
318,551,974,004 bytes. Every filename, size and SHA-256 matches across all
three workers. The retained worker completed that read in 97.95 seconds; the
replicas took 1,907.85 and 1,908.58 seconds. These checks overlap pipeline
execution and are storage/provenance work, not validation searches. Their
elapsed times must not be added to overlapping stage time as if serial.

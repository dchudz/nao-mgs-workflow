# ONT read-filtering benchmark

Compares faster variants of `FILTLONG` and `MASK_FASTQ_READS` in
`EXTRACT_VIRAL_READS_ONT`. This folder is for measurement only and isn't meant
to be merged.

## Variants

| Step | Variant | What changes |
|---|---|---|
| FILTLONG | `current` | Production command: `filtlong --verbose … \| gzip` |
| | `quiet` | Drops `--verbose`, which writes a line per read to the task log |
| | `quiet_pigz` | Also compresses output with multithreaded `pigz` |
| | `chunked` | Splits reads round-robin into one chunk per CPU and filters chunks in parallel |
| | `chopper` | Multithreaded single-pass filter with equivalent thresholds (`-q 10` ≡ `--min_mean_q 90`) |
| MASK_FASTQ_READS | `current` | Production command: `zcat \| bbmask.sh \| gzip` |
| | `pigz` | `pigz` for decompression and compression |
| | `native_gz` | BBMask reads and writes gzip itself |
| | `chunked` | Splits reads round-robin across 8 BBMask processes, each holding 1/8 of the input in memory |
| | `bbduk` | Streaming entropy masking with BBDuk instead of BBMask. Masks the same bases but sets their quality to 0, so its output differs from `current` |

Every task requests the same CPUs and memory as production, and each variant
runs `--replicates` times (default 2). Outputs of one replicate are checksummed
so each variant can be checked against `current`.

## Inputs

A samplesheet with one row per test sample:

```
label,filtlong_input,mask_input
large,<raw reads of the sample, as in the run samplesheet>,<that sample's FILTLONG output>
mid,<...>,<...>
```

`label` is a free-form name that appears in the outputs, so it shouldn't
identify the sample. `mask_input` is the `*_filtlong.fastq.gz` file from the
`FILTLONG` task's work directory of a previous run.

## Running on AWS Batch

```
nextflow run benchmarks/ont-read-filtering/main.nf -profile batch \
    --samplesheet samplesheet.csv \
    --outdir s3://<bucket>/<prefix>/results \
    --queue <on-demand job queue> \
    -w s3://<bucket>/<prefix>/work
```

Add `--batch_job_role <arn>` if the queue needs an explicit job role. Use
`--filtlong_variants` and `--mask_variants` (comma-separated) to run a subset.

## Outputs

`--outdir` gets three summary files, none of which contain reads:

- `trace.tsv`: Nextflow trace for every task. Per `.claude/benchmarking.md`,
  runtime is `complete − start` and cpu-hours is `realtime × cpus / 3600`.
- `hosts.tsv`: wall time of each task's core command, plus the vCPUs, memory and
  CPU model of the instance it ran on.
- `checksums.tsv`: read, base and masked-base (`N`) counts per variant, an
  ordered MD5 of the output, and order-independent hashes of whole records,
  read IDs, sequences, and quality scores outside `N` bases. A variant
  reproduces `current` if `hash_unordered` matches. If only
  `hash_seqs_unordered` and `hash_quals_outside_n_unordered` match, the variant
  masked the same bases and differs only in the quality scores of masked bases.

To recompute checksums after changing `bin/fastq_checksum.py`, rerun with
`-resume` and the same work directory: only the `CHECKSUM` tasks rerun.

The work directory holds full-size read files from every variant. Delete it
after the run.

## Local plumbing test

```
python3 scripts/make_synthetic_reads.py 20000 1 | gzip > synth.fastq.gz
# samplesheet with synth.fastq.gz in both input columns, then:
nextflow run main.nf -profile local --samplesheet samplesheet.csv \
    --outdir results --replicates 1
```

The synthetic reads are random, so they test correctness and plumbing but not
speed.

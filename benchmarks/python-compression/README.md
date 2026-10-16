# Python output-compression benchmark

Measures how much of the slowest Illumina Python steps goes to compressing
their output. Every module script writes `.gz` files through Python's
`gzip.open`, which defaults to compression level 9. This folder is for
measurement only and isn't meant to be merged.

## How it works

Each step runs a production script from `modules/local/` unchanged, with the
same container, CPUs and memory as production. `hook/sitecustomize.py` is put on
`PYTHONPATH`, so Python loads it at startup. It replaces `gzip.open` for writes
with the writer that `BENCH_COMPRESSION` names. Reads are untouched, so every
variant reads the same production gzip inputs.

| Variant | Output written as |
|---|---|
| `gzip9` | gzip level 9, the `gzip.open` default (production) |
| `gzip6` | gzip level 6, the `gzip` command-line default |
| `gzip1` | gzip level 1 |
| `zstd3` | Zstandard level 3, the `zstd` default. Needs Python 3.14, so it's skipped for `PROCESS_VIRAL_BOWTIE2_SAM`, whose container has 3.13 |
| `none` | uncompressed: the floor, showing the step's time without compression |

Outputs keep their `.gz` names whatever the format. Each variant runs
`--replicates` times (default 2).

## Steps and inputs

Each step takes the inputs of one task from a run, ideally the slowest one in
a large delivery. Find the task's work directory in the run's trace. Its
`.command.sh` shows the command with the input names in argument order, and its
`.command.run` shows where each input was staged from. If the run's work files
are gone, rerun `RUN` on just that library (one samplesheet row, the run's
params) and take the inputs from the new run.

| Step | Production process | Params | Input |
|---|---|---|---|
| `FILTER_VIRAL_SAM` | same | `fvs_sam`, `fvs_reads` | The sorted SAM (first argument) and the sorted reads (second) |
| `PROCESS_VIRAL_BOWTIE2_SAM` | same | `pvbs_sam` | The `FILTER_VIRAL_SAM` output, `*_viral_filtered.sam.gz` |
| `ADD_SAMPLE_COLUMN` | `ADD_SAMPLE_COLUMN_ALIGNER` | `asc_tsv` | The `PROCESS_VIRAL_BOWTIE2_SAM` output |
| `JOIN_TSVS` | `PROCESS_LCA_ALIGNER_OUTPUT:JOIN_TSVS` | `jt_tsv1`, `jt_tsv2` | The task's two inputs, in argument order |
| `FILTER_TSV_COLUMN_BY_VALUE` | same | `ftcbv_tsv` | The `JOIN_TSVS` output |
| `CHECK_TSV_DUPLICATES` | same (downstream) | `ctd_tsv` | The task's input |
| `ADD_CONDITIONAL_TSV_COLUMN` | same (downstream) | `act_tsv` | The downstream `JOIN_TSVS` output |

The first five run in that order on one library in `RUN`, so one rerun gives
all their inputs. For downstream's `JOIN_TSVS` instead, also set
`jt_field=aligner_taxid_lca` and `jt_type=left`.

Also set `ref_dir` to the run's index (for `PROCESS_VIRAL_BOWTIE2_SAM`'s
metadata files). `bt2_score_threshold` defaults to production's 20; check it
against `FILTER_VIRAL_SAM`'s `.command.sh`. A step whose params are unset is
skipped.

Put the params in a JSON file and keep it out of the repository:

```
{
  "ref_dir": "<the run's ref_dir>",
  "fvs_sam": "<...>",
  "fvs_reads": "<...>",
  "pvbs_sam": "<...>",
  "asc_tsv": "<...>",
  "jt_tsv1": "<...>",
  "jt_tsv2": "<...>",
  "ftcbv_tsv": "<...>"
}
```

## Running on AWS Batch

```
nextflow run benchmarks/python-compression/main.nf -profile batch \
    -params-file inputs.json \
    --outdir s3://<bucket>/<prefix>/results \
    --queue <on-demand job queue> \
    -w s3://<bucket>/<prefix>/work
```

Add `--batch_job_role <arn>` if the queue needs an explicit job role. Use
`--variants` (comma-separated) to run a subset.

## Outputs

`--outdir` gets three summary files, none of which contain reads:

- `trace.tsv`: Nextflow trace for every task. Per `.claude/benchmarking.md`,
  runtime is `complete − start` and cpu-hours is `realtime × cpus / 3600`.
- `hosts.tsv`: wall time of each task's core command, plus the vCPUs, memory and
  CPU model of the instance it ran on.
- `checksums.tsv`: for one replicate of each step and variant, the output's
  detected format, compressed and decompressed size, line count and MD5 of the
  decompressed content. Within a step, every variant should match `gzip9` on
  everything but format and compressed size; the format should match the
  variant, which confirms the hook was active.

The work directory holds full-size outputs from every variant. Delete it after
the run.

## Local plumbing test

```
python3 scripts/make_toy_inputs.py toy
nextflow run main.nf -profile local -params-file toy/toy-params.json \
    --outdir results --replicates 1
```

The toy inputs are synthetic and small, so they test correctness and plumbing
but not speed.

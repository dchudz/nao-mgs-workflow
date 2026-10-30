Written by Codex, an AI, at David's request; David did not review it line by line.

# Experiment sources and archived outputs

This branch retains experiment code, configurations, synthetic test fixtures, and report figures. Generated outputs, execution logs, result tables, and superseded report drafts are omitted from the final diff. No experimental code or results were recomputed by this cleanup.

The complete pre-cleanup tree remains at [commit 7035f902](https://github.com/dchudz/nao-mgs-workflow/tree/7035f90215403171fe82425860668d776dfad0b2/experiments). The cleanup adds a commit; it does not rewrite experiment history. Every removed file was verified against its published Git blob before removal.

The [PR 5 description](https://github.com/dchudz/nao-mgs-workflow/pull/5) contains the current report, a guide to key experiment commits, and private S3 locations with preservation receipts. Historical README files describe their original experiment state and may refer to outputs that are present only in the archive.

For analysis scripts or figure generators that read archived output paths, extract the original tree into a separate directory first:

```sh
mkdir -p /tmp/pr5-archived-experiments
git archive 7035f90215403171fe82425860668d776dfad0b2 experiments | tar -x -C /tmp/pr5-archived-experiments
```

Synthetic test inputs and expected answers under `paired_pipeline_comparison/synthetic-control-fixture/` remain tracked. They are test fixtures, not the real-data result tables.

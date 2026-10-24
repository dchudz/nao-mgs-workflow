# Local archive-compression control

One existing checksum-verified archive was decompressed to a fixed 1,507,624,960-byte tar payload. Each compression level was run twice in reversed order; all six outputs restored exactly the same payload hash. This exercised preservation only. No sequences were interpreted and no running pipeline or existing backup was changed.

| gzip level | Median compression and fsync time | Compressed bytes | Time relative to level 9 | Size relative to level 9 |
|---|---:|---:|---:|---:|
| 1 | 3.666 s | 99,599,124 | 0.147 | 1.411 |
| 6 | 11.930 s | 77,400,331 | 0.477 | 1.097 |
| 9 | 24.990 s | 70,574,238 | 1.000 | 1.000 |

These are local-computer observations on one archive, not measured AWS or end-to-end speedups. The input was already cached. Timing excludes the subsequent decompression check. All benchmark gzip headers used an empty filename and zero timestamp, so even level 9 has a different outer checksum from the original archive. The restored tar bytes match exactly.

The executed source and full measurements are saved beside this file. To reproduce in a scratch directory, copy `executed-script.py` there and place the source archive and its receipt at `interrupted-worker-completed-evidence.tar.gz` and `interrupted-worker-completed-backup.json`. The source S3 URI, compressed checksum, and restored-payload checksum are recorded in `result.json`. Running the script creates its outputs under `compression-control/`.

The running cohort continues to use its frozen archive implementation. These measurements support only the stated local compression/size tradeoff.

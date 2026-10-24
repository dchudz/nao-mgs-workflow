# Local lossless archive codec measurement

This records a completed local preservation experiment, not a pipeline change.
The same 1,507,624,960-byte tar payload used in the gzip control was compressed
with Python 3.15 and Zstandard 1.5.7 at level 3. Both outputs were decompressed
and their full payload SHA-256 matched the original.

| Codec | Median compression plus fsync seconds | Output bytes |
| --- | ---: | ---: |
| gzip level 9 | 24.990 | 70,574,238 |
| gzip level 6 | 11.930 | 77,400,331 |
| gzip level 1 | 3.666 | 99,599,124 |
| Zstandard level 3 | 1.368 | 79,477,184 |

There were two observations per setting on one cached local payload. These
measurements do not establish AWS performance, transfer savings, or end-to-end
pipeline speed. Zstandard requires a compatible reader; deployment compatibility
was not tested. The running study and original archives were unchanged.

`result.json` contains exact measurements and checksums. `executed-script.py`
is the byte-identical executed script. To reproduce, place it in a scratch
directory containing `compression-control/result.json` and the verified
`compression-control/payload.tar` used by the preceding experiment. Do not run
it directly in this results directory. The preceding gzip result is retained
in the sibling `compression-control` directory.

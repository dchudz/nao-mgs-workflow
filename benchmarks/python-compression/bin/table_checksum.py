#!/usr/bin/env python3
"""Summarize one benchmark output file, whatever its compression.

Detects gzip, Zstandard or no compression from the file's first bytes, then
streams the decompressed content through MD5. Prints one TSV row (with header):
step, variant, format, compressed bytes, decompressed bytes, line count and MD5.
Variants of a step should agree on everything except format and compressed size.
"""

import argparse
import gzip
import hashlib
import os
from typing import IO

GZIP_MAGIC = b"\x1f\x8b"
ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
CHUNK = 1 << 20


def detect_format(path: str) -> str:
    """Return "gzip", "zstd" or "none" from the file's magic bytes."""
    with open(path, "rb") as f:
        head = f.read(4)
    if head.startswith(GZIP_MAGIC):
        return "gzip"
    if head.startswith(ZSTD_MAGIC):
        return "zstd"
    return "none"


def open_decompressed(path: str, fmt: str) -> IO[bytes]:
    """Open path for reading decompressed bytes."""
    if fmt == "gzip":
        return gzip.open(path, "rb")
    if fmt == "zstd":
        from compression import zstd  # Python 3.14+

        return zstd.open(path, "rb")
    return open(path, "rb")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step")
    parser.add_argument("variant")
    parser.add_argument("path")
    args = parser.parse_args()

    fmt = detect_format(args.path)
    md5 = hashlib.md5()
    n_bytes = 0
    n_lines = 0
    with open_decompressed(args.path, fmt) as f:
        while chunk := f.read(CHUNK):
            md5.update(chunk)
            n_bytes += len(chunk)
            n_lines += chunk.count(b"\n")

    header = ["step", "variant", "format", "compressed_bytes", "bytes", "lines", "md5"]
    row = [
        args.step,
        args.variant,
        fmt,
        os.path.getsize(args.path),
        n_bytes,
        n_lines,
        md5.hexdigest(),
    ]
    print("\t".join(header))
    print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()

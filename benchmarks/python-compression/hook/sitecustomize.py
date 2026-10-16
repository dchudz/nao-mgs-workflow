"""Benchmark hook: change how Python scripts compress the .gz files they write.

Python imports sitecustomize at startup from any directory on PYTHONPATH, so the
production scripts run unchanged. BENCH_COMPRESSION picks the writer that
gzip.open() uses in write, append or exclusive-create mode:

    gzip9  gzip at level 9, the gzip.open default (production)
    gzipN  gzip at level N (1-9)
    zstdN  Zstandard at level N (needs Python 3.14+)
    none   no compression

Files keep their .gz names whatever the format, and reads are untouched.
Without BENCH_COMPRESSION the hook does nothing.

Python reports and then ignores errors raised while importing sitecustomize, so
the variant is parsed at the first write instead: a bad variant fails the script.
"""

import builtins
import gzip
import os
import re
import sys
from typing import IO, Any

_VARIANT = os.environ.get("BENCH_COMPRESSION")
_gzip_open = gzip.open


def _open_for_write(
    filename: Any,
    mode: str,
    encoding: str | None,
    errors: str | None,
    newline: str | None,
) -> IO[Any]:
    """Open filename for writing with the compression BENCH_COMPRESSION names."""
    assert _VARIANT is not None
    if _VARIANT == "none":
        if "t" not in mode and "b" not in mode:
            mode += "b"  # gzip.open treats a bare "w" as binary
        return builtins.open(
            filename, mode, encoding=encoding, errors=errors, newline=newline
        )
    match = re.fullmatch(r"(gzip|zstd)([0-9]+)", _VARIANT)
    if not match:
        raise ValueError(f"Unknown BENCH_COMPRESSION: {_VARIANT!r}")
    level = int(match.group(2))
    if match.group(1) == "gzip":
        return _gzip_open(filename, mode, level, encoding, errors, newline)
    from compression import zstd  # Python 3.14+

    return zstd.open(
        filename, mode, level=level, encoding=encoding, errors=errors, newline=newline
    )


def _bench_open(
    filename: Any,
    mode: str = "rb",
    compresslevel: int = 9,
    encoding: str | None = None,
    errors: str | None = None,
    newline: str | None = None,
) -> IO[Any]:
    """Stand-in for gzip.open that redirects writes to the benchmark variant."""
    if "r" in mode:
        return _gzip_open(filename, mode, compresslevel, encoding, errors, newline)
    print(f"[bench-compression] {_VARIANT}: {filename}", file=sys.stderr)
    return _open_for_write(filename, mode, encoding, errors, newline)


if _VARIANT:
    gzip.open = _bench_open  # type: ignore[assignment]

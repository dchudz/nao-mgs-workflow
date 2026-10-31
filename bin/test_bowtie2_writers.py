# Written by Codex, an AI, at David's request; David has not reviewed it line by line.
"""Check the BOWTIE2 shell lifecycle without downloads, alignments, or Nextflow.

Only fixed fixture parameters are rendered; this does not validate Groovy syntax.
Command doubles pass inert bytes through the real Bash/tee/FIFO topology. The
compressor double uses Python gzip and can fail after consuming all its input.
Run with: python3 -m unittest discover -s bin -p test_bowtie2_writers.py
"""

import gzip
import os
import re
import signal
import subprocess
import sys
import tempfile
import unittest
from contextlib import suppress
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = b"writer lifecycle fixture\n"
STUB = r"""#!/usr/bin/env python3
import gzip
import os
from pathlib import Path
import sys
import time

name = Path(sys.argv[0]).name
args = sys.argv[1:]
if name == "download_db.py":
    print(".")
elif name == "pigz" and "-dc" in args:
    sys.stdout.buffer.write(gzip.decompress(Path(args[-1]).read_bytes()))
else:
    target = Path(os.environ["TARGET"])
    selected = (name in ("pigz", "gzip") and target.exists()
                and os.path.samestat(target.stat(), os.fstat(1)))
    scenario = os.environ["SCENARIO"]
    if selected and scenario == "early_failure":
        sys.exit(23)
    data = sys.stdin.buffer.read()
    if name == "bowtie2" and scenario == "producer_failure":
        sys.exit(17)
    if name == "samtools" and "-f" in args and scenario == "intermediate_failure":
        sys.exit(19)
    if selected:
        if scenario in ("delay", "producer_failure"):
            time.sleep(0.5)
        if scenario == "writer_failure":
            sys.exit(23)
    sys.stdout.buffer.write(gzip.compress(data, mtime=0)
                            if name in ("pigz", "gzip") else data)
"""


def render(debug: bool) -> str:
    source = (ROOT / "modules/local/bowtie2/main.nf").read_text()
    values = {
        "index_dir": ".",
        "params_map.db_download_timeout": "1",
        "task.cpus": "1",
        "par": "",
        "unmapped_flag": "12",
        "reads_interleaved": "reads.gz",
        "sam": "mapped.sam.gz",
        "al": "mapped.fastq.gz",
        "un": "unmapped.fastq.gz",
        "in2": "input.fastq.gz",
        "params_map.debug": str(debug).lower(),
    }
    # Read the debug fragments from the module rather than copying their values.
    for name, fragment in re.findall(
        r"def (debug_\w+) = params_map.debug \? '([^']*)' : \"\"", source
    ):
        values[name] = fragment if debug else ""
    body = source.split('"""')[1].replace(r"\$", "__DOLLAR__")
    for name, value in values.items():
        body = body.replace("${" + name + "}", value.replace("$", "__DOLLAR__"))
    ternary = re.compile(
        r"\$\{\s*params_map\.(debug|interleaved|remove_sq)\s*\?\s*"
        r'"((?:\\.|[^"\\])*)"\s*:\s*"((?:\\.|[^"\\])*)"\s*\}'
    )
    body = ternary.sub(
        lambda m: (
            m[2]
            if {"debug": debug, "interleaved": True, "remove_sq": False}[m[1]]
            else m[3]
        ),
        body,
    )
    if "${" in body:
        raise ValueError("Unsupported fixture expression in module shell")
    return re.sub(r'\\([\\"])', r"\1", body).replace("__DOLLAR__", "$")


class WriterLifecycle(unittest.TestCase):
    def test_all_writers_finish_and_failures_reach_the_caller(self) -> None:
        cases = [
            ("success", "unmapped.fastq.gz"),
            ("writer_failure", "unmapped.fastq.gz"),
            ("writer_failure", "mapped.fastq.gz"),
            ("writer_failure", "mapped.sam.gz"),
            ("early_failure", "unmapped.fastq.gz"),
            ("output_open_failure", "unmapped.fastq.gz"),
            ("intermediate_failure", "unmapped.fastq.gz"),
            ("producer_failure", "unmapped.fastq.gz"),
        ]
        if Path("/dev/full").exists():
            cases.append(("disk_full", "unmapped.fastq.gz"))
        for debug in (False, True):
            debug_cases = (
                [
                    (scenario, name)
                    for name in (
                        "test_all.sam.gz",
                        "test_unmapped.sam.gz",
                        "test_mapped.sam.gz",
                    )
                    for scenario in ("writer_failure", "delay")
                ]
                if debug
                else []
            )
            for scenario, target in cases + debug_cases:
                with self.subTest(debug=debug, scenario=scenario, target=target):
                    self.check_case(debug, scenario, target)

    def check_case(self, debug: bool, scenario: str, target: str) -> None:
        with tempfile.TemporaryDirectory(prefix="bowtie-writers-") as directory:
            wd = Path(directory)
            (wd / "bin").mkdir()
            (wd / "tmp").mkdir()
            stub = wd / "bin" / "stub"
            stub.write_text(STUB.replace("/usr/bin/env python3", sys.executable))
            stub.chmod(0o755)
            for name in (
                "download_db.py",
                "pigz",
                "gzip",
                "bowtie2",
                "samtools",
                "sed",
            ):
                (wd / "bin" / name).symlink_to("stub")
            payload = PAYLOAD * (8192 if scenario == "early_failure" else 1)
            (wd / "reads.gz").write_bytes(gzip.compress(payload))
            if scenario == "output_open_failure":
                (wd / target).mkdir()
            elif scenario == "disk_full":
                (wd / target).symlink_to("/dev/full")
            env = dict(
                os.environ,
                PATH=str(wd / "bin") + ":" + os.environ["PATH"],
                TMPDIR=str(wd / "tmp"),
                SCENARIO=scenario,
                TARGET=target,
            )
            (wd / "command.sh").write_text(render(debug))
            with (wd / "stderr").open("wb") as err:
                proc = subprocess.Popen(
                    ["bash", "command.sh"],
                    cwd=wd,
                    env=env,
                    start_new_session=True,
                    stdout=subprocess.DEVNULL,
                    stderr=err,
                )
                try:
                    code = proc.wait(timeout=10)
                    diagnostic = (wd / "stderr").read_text()
                    if scenario in ("success", "delay"):
                        self.assertEqual(code, 0, diagnostic)
                        names = [
                            "mapped.sam.gz",
                            "mapped.fastq.gz",
                            "unmapped.fastq.gz",
                        ]
                        if debug:
                            names += [
                                "test_all.sam.gz",
                                "test_unmapped.sam.gz",
                                "test_mapped.sam.gz",
                            ]
                        for name in names:
                            raw = (wd / name).read_bytes()
                            self.assertTrue(raw, name + " was empty at task exit")
                            self.assertEqual(gzip.decompress(raw), payload, name)
                    else:
                        self.assertNotEqual(
                            code, 0, "Writer failure was reported as success"
                        )
                    with self.assertRaises(
                        ProcessLookupError, msg="Child outlived task"
                    ):
                        os.killpg(proc.pid, 0)
                    self.assertEqual(
                        list((wd / "tmp").iterdir()), [], "FIFO scratch leaked"
                    )
                finally:
                    with suppress(ProcessLookupError):
                        os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()


if __name__ == "__main__":
    unittest.main()

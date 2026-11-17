'''
Run the mapper on every example set and print one summary table.

python run_tests.py                     # all sets, mapper = main.py
python run_tests.py --quick             # skip the 20 Mbp sets
python run_tests.py --mapper mapper.py  # test another entry point

Besides the given files it builds a "mixed" set (reads1 + reads_unmappable, the unmappable
reads renamed so the ids don't clash), and checks that output lines are tab-separated.
Exit code is 1 if any minimum requirement fails.
'''

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from evaluate import compare, count_reads, pct, read_coords

ROOT = Path(__file__).resolve().parent
D1 = ROOT / "data" / "example_files1"
D2 = ROOT / "data" / "example_files2"


def make_mixed(tmp):
    """reads1 + reads_unmappable (ids prefixed with 'unmappable_'); truth = reads1.txt."""

    mixed = tmp / "reads_mixed.fasta"
    with open(mixed, "w") as out:
        out.write((D1 / "reads1.fasta").read_text())
        for line in (D1 / "reads_unmappable.fasta").read_text().splitlines(keepends=True):
            if line.startswith(">"):
                line = ">unmappable_" + line[1:].strip() + "\n"
            out.write(line)
    return mixed


def run_mapper(mapper, ref, reads, output):
    """Run the mapper, return (elapsed seconds, peak memory MB, exit code, stderr)."""

    err_path = output.with_suffix(".stderr")
    with open(err_path, "w") as err:
        t0 = time.perf_counter()
        proc = subprocess.Popen([sys.executable, mapper, ref, reads, output],
                                stdout=subprocess.DEVNULL, stderr=err)
        # wait4 gives the memory of this one child (RUSAGE_CHILDREN would be the max of all runs)
        _, status, usage = os.wait4(proc.pid, 0)
        elapsed = time.perf_counter() - t0

    # ru_maxrss is in bytes on macOS and in kilobytes on Linux
    memory = usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    return elapsed, memory, os.waitstatus_to_exitcode(status), err_path.read_text()


def tabs_ok(path):
    """Every non-empty output line must be exactly id<TAB>start<TAB>end."""

    with open(path) as f:
        return all(line.rstrip("\n").count("\t") == 2 for line in f if line.strip())


def main():
    parser = argparse.ArgumentParser(description="Run the mapper on all example sets")
    parser.add_argument("--mapper", default=str(ROOT / "main.py"), help="mapper script (default main.py)")
    parser.add_argument("--quick", action="store_true", help="skip the 20 Mbp sets")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # (name, reference, reads, truth or None for "nothing should map")
        cases = [
            ("reads0", D1 / "reference.fasta", D1 / "reads0.fasta", D1 / "reads0.txt"),
            ("reads1", D1 / "reference.fasta", D1 / "reads1.fasta", D1 / "reads1.txt"),
            ("reads2", D1 / "reference.fasta", D1 / "reads2.fasta", D1 / "reads2.txt"),
            ("unmappable", D1 / "reference.fasta", D1 / "reads_unmappable.fasta", None),
            ("mixed", D1 / "reference.fasta", make_mixed(tmp), D1 / "reads1.txt"),
        ]
        if not args.quick:
            cases += [
                ("reads20Ma", D2 / "reference20M.fasta", D2 / "reads20Ma.fasta", D2 / "reads20Ma.txt"),
                ("reads20Mb", D2 / "reference20M.fasta", D2 / "reads20Mb.fasta", D2 / "reads20Mb.txt"),
            ]

        header = f"{'set':<11}{'reads':>6}{'mappable':>9}{'mapped%':>9}{'wrong':>6}{'wrong%':>8}" \
                 f"{'max dev':>8}{'s/read':>8}{'total s':>8}{'MB':>6}{'tabs':>6}  result"
        print(header)
        print("-" * len(header))

        all_ok = True
        for name, ref, reads, truth in cases:
            output = tmp / f"out_{name}.txt"
            elapsed, memory, code, stderr = run_mapper(args.mapper, ref, reads, output)
            if code != 0:
                print(f"{name:<11} mapper failed (exit code {code}):\n{stderr}")
                all_ok = False
                continue

            n_reads = count_reads(reads)
            out = read_coords(output)
            expected = read_coords(truth) if truth else {}
            correct, wrong, mapped_mappable = compare(out, expected)
            mapped_pct = pct(mapped_mappable, len(expected))
            wrong_pct = pct(len(wrong), len(out))
            tabs = tabs_ok(output)

            failed = []
            if expected and mapped_pct < 80:
                failed.append("mapped<80%")
            if wrong_pct > 1:
                failed.append("wrong>1%")
            if elapsed > 100 + 3 * n_reads:
                failed.append("time")
            if memory >= 1024:
                failed.append("memory")
            if not tabs:
                failed.append("tabs")
            all_ok = all_ok and not failed

            print(f"{name:<11}{n_reads:>6}{len(expected):>9}"
                  f"{(f'{mapped_pct:.1f}' if expected else '-'):>9}{len(wrong):>6}{wrong_pct:>8.2f}"
                  f"{(max(correct) if correct else '-'):>8}{elapsed / n_reads:>8.3f}{elapsed:>8.1f}"
                  f"{memory:>6.0f}{('ok' if tabs else 'BAD'):>6}  {'FAIL: ' + ', '.join(failed) if failed else 'OK'}")
            for read_id, reason in wrong[:5]:
                print(f"{'':<11}wrong: {read_id} - {reason}")

    print("\nall minimum requirements met" if all_ok else "\nsome minimum requirements are NOT met")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()

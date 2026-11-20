"""
How to run the code:

python evaluate.py reads1.fasta output.txt reads1.txt
python evaluate.py reads_unmappable.fasta output.txt
python evaluate.py reads1.fasta output.txt reads1.txt --mapper main.py --ref reference.fasta
"""


import argparse
import subprocess
import sys
import time

tolerance = 20 


def count_reads(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def read_coords(path):
    coords = {}
    with open(path) as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if line == "":
                continue
            parts = line.split()
            if len(parts) != 3:
                print(f"WARNING: {path}, line {line_no}: expected 3 fields, got {len(parts)} - skipped")
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                print(f"WARNING: {path}, line {line_no}: start/end is not a number - skipped")
                continue
            if parts[0] in coords:
                print(f"WARNING: read {parts[0]} appears more than once in {path}")
            coords[parts[0]] = (start, end)
    return coords


def run_mapper(mapper, ref, reads, output):
    cmd = [sys.executable, mapper, ref, reads, output]
    print("running:", " ".join(cmd))
    t0 = time.perf_counter()
    result = subprocess.run(cmd)
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print("the mapper failed, exit code", result.returncode)
        sys.exit(1)
    return elapsed, peak_memory_mb()


def peak_memory_mb():
    try:
        import resource
    except ImportError:
        return None
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == "darwin":
        return peak / 1024 / 1024
    return peak / 1024


def compare(output, expected):
    correct = []
    wrong = []
    mapped_mappable = 0
    for read_id in output:
        start, end = output[read_id]
        if read_id not in expected:
            wrong.append((read_id, "this read should not be mapped"))
            continue
        mapped_mappable += 1
        exp_start, exp_end = expected[read_id]
        d_start = abs(start - exp_start)
        d_end = abs(end - exp_end)
        if d_start <= tolerance and d_end <= tolerance:
            correct.append(max(d_start, d_end))
        else:
            wrong.append((read_id, f"got {start}-{end}, expected {exp_start}-{exp_end}"))
    return correct, wrong, mapped_mappable


def pct(a, b):
    if b == 0:
        return 0.0
    return 100.0 * a / b


def points_correctness(wrong_percent):
    if wrong_percent == 0:
        return 3
    if wrong_percent <= 0.2:
        return 2
    if wrong_percent <= 0.5:
        return 1
    return 0


def points_completeness(mapped_percent):
    if mapped_percent >= 99:
        return 3
    if mapped_percent >= 95:
        return 2
    if mapped_percent >= 90:
        return 1
    return 0


def points_time(sec_per_read):
    if sec_per_read <= 0.2:
        return 3
    if sec_per_read <= 0.5:
        return 2
    if sec_per_read <= 1:
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description="Evaluate the output of the read mapper")
    parser.add_argument("reads", help="reads fasta file")
    parser.add_argument("output", help="mapper output file (created if --mapper is used)")
    parser.add_argument("expected", nargs="?", help="file with true coordinates (leave out for unmappable reads)")
    parser.add_argument("--mapper", help="path to main.py, it will be run before the evaluation")
    parser.add_argument("--ref", help="reference fasta (needed with --mapper)")
    parser.add_argument("--time", type=float, help="running time in seconds (if you ran the mapper by hand)")
    args = parser.parse_args()

    elapsed = args.time
    memory = None
    if args.mapper:
        if not args.ref:
            parser.error("--mapper needs --ref")
        elapsed, memory = run_mapper(args.mapper, args.ref, args.reads, args.output)

    n_reads = count_reads(args.reads)
    output = read_coords(args.output)
    expected = {}
    if args.expected:
        expected = read_coords(args.expected)

    correct, wrong, mapped_mappable = compare(output, expected)
    n_out = len(output)
    n_mappable = len(expected)

    print()
    print("reads in the fasta file:   ", n_reads)
    print("mappable reads (expected): ", n_mappable)
    print("reads in the output:       ", n_out)
    if n_mappable > 0:
        print(f"mapped (of mappable): {mapped_mappable} ({pct(mapped_mappable, n_mappable):.2f}%)")
        print(f"mapped correctly: {len(correct)} ({pct(len(correct), n_mappable):.2f}% of mappable)")
    print(f"mapped incorrectly: {len(wrong)} ({pct(len(wrong), n_out):.2f}% of output lines, "
          f"{pct(len(wrong), n_reads):.2f}% of all reads)")
    if len(correct) > 0:
        print(f"deviation of correct reads:  max {max(correct)} bp, mean {sum(correct) / len(correct):.1f} bp")
 
    if len(wrong) > 0:
        print("\nfirst wrong reads:")
        for read_id, reason in wrong[:10]:
            print("  ", read_id, "-", reason)
        if len(wrong) > 10:
            print("   ...", len(wrong) - 10, "more")
 
    sec_per_read = None
    if elapsed is not None and n_reads > 0:
        sec_per_read = elapsed / n_reads
        limit = 100 + 3 * n_reads
        print(f"\ntime: {elapsed:.1f} s in total, {sec_per_read:.3f} s per read (limit {limit} s)")
    if memory is not None:
        print(f"peak memory of the mapper: {memory:.0f} MB (limit 1024 MB)")
 
    print("\nminimum requirements:")
    ok_all = True
    if n_mappable > 0:
        ok = pct(mapped_mappable, n_mappable) >= 80
        ok_all = ok_all and ok
        print("  >= 80% of mappable reads mapped:", "OK" if ok else "FAILED")
    ok = pct(len(wrong), n_out) <= 1
    ok_all = ok_all and ok
    print("  <= 1% mapped incorrectly:        ", "OK" if ok else "FAILED")
    if elapsed is not None:
        ok = elapsed <= 100 + 3 * n_reads
        ok_all = ok_all and ok
        print("  time <= 100 + 3r seconds:        ", "OK" if ok else "FAILED")
    else:
        print("  time: unknown (use --time or --mapper)")
    if memory is not None:
        ok = memory < 1024
        ok_all = ok_all and ok
        print("  memory < 1 GB:                   ", "OK" if ok else "FAILED")
    else:
        print("  memory: unknown (measured only with --mapper, on Linux)")
 
    print("\nestimated points:")
    p = points_correctness(pct(len(wrong), n_out))
    print("  correctness: ", p, "/ 3")
    if n_mappable > 0:
        print("  completeness:", points_completeness(pct(mapped_mappable, n_mappable)), "/ 3")
    if sec_per_read is not None:
        print("  time:        ", points_time(sec_per_read), "/ 3")
    if not ok_all:
        print("\nsome of the minimum requirements are not met!")
 
 
if __name__ == "__main__":
    main()

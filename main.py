"""
Read mapper — main file
"""

import argparse
import logging
import os
import time

from modules.config import K
from modules.parse_encode import read_reference, read_fasta, encode_str
from modules.index import build_index
from modules.decision import map_read


def setup_logging(output_path):
    """Log to stdout and to <output>.log."""

    log_path = os.path.splitext(output_path)[0] + ".log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_path, mode="w")],
    )


def main():

    parser = argparse.ArgumentParser(description="Seed-and-extend read mapper")
    parser.add_argument("--reference", required=True, help="Reference FASTA file (<= 20 Mbp)")
    parser.add_argument("--reads", required=True, help="FASTA file with reads you want to map (~1 kbp each)")
    parser.add_argument("--output", required=True, help="Output file, one line per mapped read")
    args = parser.parse_args()

    setup_logging(args.output)
    t_start = time.perf_counter()

    # Index the reference file.
    
    ref_encoded = encode_str(read_reference(args.reference))
    logging.info("Reference loaded: %d bp", len(ref_encoded))
    t0 = time.perf_counter()

    index = build_index(ref_encoded, K)
    logging.info("Index built in %.2f s", time.perf_counter() - t0)

    # Steps 4–9: map each read
    n_reads = n_mapped = 0
    t0 = time.perf_counter()
    with open(args.output, "w") as out:
        for read_id, seq in read_fasta(args.reads):
            n_reads += 1
            result = map_read(encode(seq), index, ref_enc)
            if result is not None:
                start, end = result
                out.write(f"{read_id}\t{start}\t{end}\n")
                n_mapped += 1
    t_map = time.perf_counter() - t0

    logging.info("Mapped %d / %d reads", n_mapped, n_reads)
    if n_reads:
        logging.info("Mapping time: %.2f s (%.4f s/read)", t_map, t_map / n_reads)
    logging.info("Total time: %.2f s", time.perf_counter() - t_start)


if __name__ == "__main__":
    main()

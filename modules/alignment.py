"""Module for finding hits and alignment of the reads."""

import numpy as np
from modules.config import MARGIN, MAX_OCC, K
from modules.index import hash_kmers


def find_hits(read_enc: np.array, index: tuple[np.array, np.array]) -> tuple[np.array, np.array]:
    """Return arrays (ref_pos, read_pos) of all exact seed hits.

    k-mers with more than MAX_OCC occurrences in the reference are ignored.
    """

    index_codes, positions = index
    read_codes, valid = hash_kmers(read_enc, K)

    # Remember where each k-mer starts in the read before dropping the invalid ones
    read_kmer_pos = np.nonzero(valid)[0]
    read_codes = read_codes[valid]

    # All occurrences of read_codes[i] in the reference are positions[lo[i]:hi[i]]
    lo = np.searchsorted(index_codes, read_codes, side="left")
    hi = np.searchsorted(index_codes, read_codes, side="right")
    counts = hi - lo # number of occurences of i k-mer.

    # Filter k-mers occuring too often in the reference, not really informative (MAX_OCC) or not occuring at all.
    keep = (counts > 0) & (counts <= MAX_OCC)
    lo, counts, read_kmer_pos = lo[keep], counts[keep], read_kmer_pos[keep]

    read_pos, idx = [], []
    for p, l, c in zip(read_kmer_pos, lo, counts):
        # If we have multiple k-mer hits in the reference make a list where every hit is represented in it
        for t in range(c):
            read_pos.append(p)
            idx.append(l + t)
    ref_pos = positions[idx]

    return np.array(ref_pos), np.array(read_pos)

def semiglobal_dp(read, window):
    """Fill the DP matrix (zero first row); return (D, edit_distance, end_col)."""
    raise NotImplementedError


def traceback(D, read, window, end_col):
    """Return the start column of the alignment ending at end_col."""
    raise NotImplementedError


def verify(candidate, ref_enc):
    """Build the reference window, run DP, return (edit_distance, start, end).

    candidate is (read_enc, est_start); start/end are reference coordinates.
    """
    raise NotImplementedError


def best_diagonals(ref_pos, read_pos):
    """Return the top candidate windows as a list of (est_start, votes).

    The first entry is the densest BAND-wide window of diagonals, the second
    (if any) is the best non-overlapping runner-up.
    """
    raise NotImplementedError

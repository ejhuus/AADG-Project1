"""Module for finding hits and alignment of the reads."""

import numpy as np
from modules.config import BAND, MARGIN, MAX_OCC, K
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

    # Expand every k-mer into one hit per reference occurrence, i.e. k-mer i contributes
    # index entries lo[i], lo[i]+1, ..., lo[i]+counts[i]-1, all paired with read_kmer_pos[i].
    # offsets[i] = where k-mer i's block starts in the flat output (exclusive cumulative sum).
    offsets = np.cumsum(counts) - counts
    read_pos = np.repeat(read_kmer_pos, counts)
    # Flat output slot h belongs to k-mer i: idx = lo[i] + (h - offsets[i])
    idx = np.repeat(lo - offsets, counts) + np.arange(counts.sum())
    ref_pos = positions[idx]

    return ref_pos, read_pos

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


def best_diagonals(ref_pos: np.array, read_pos: np.array) -> list[tuple[int, int]]:
    """Return the top candidate windows as a list of (est_start, votes).

    The first entry is the densest BAND-wide window of diagonals, the second
    (if any) is the best non-overlapping runner-up.
    """

    if len(ref_pos) == 0:
        return []

    # Each hit votes for the diagonal ref_pos - read_pos (= estimated read start in the reference).
    # int64 because the positions are unsigned and the difference can be negative.
    diag = np.sort(ref_pos.astype(np.int64) - read_pos.astype(np.int64))

    # Window i covers diagonals [diag[i], diag[i] + BAND); ends[i] is the first hit past it,
    # so votes[i] = number of hits in window i.
    ends = np.searchsorted(diag, diag + BAND, side="left")
    votes = ends - np.arange(len(diag))

    def candidate(i):
        # Median diagonal of the window: robust to the indel drift at the edges.
        return int(diag[(i + ends[i] - 1) // 2]), int(votes[i])

    best = int(np.argmax(votes))
    candidates = [candidate(best)]

    # Runner-up: best window that does not overlap the best one.
    lo, hi = diag[best], diag[best] + BAND
    outside = (diag + BAND <= lo) | (diag >= hi)
    if outside.any():
        second = int(np.argmax(np.where(outside, votes, 0)))
        candidates.append(candidate(second))

    return candidates

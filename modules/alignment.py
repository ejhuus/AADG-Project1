"""Module for finding hits and alignment of the reads."""

import numpy as np

from modules.config import BAND, MARGIN, MAX_ERR_RATE, MAX_OCC, MIN_ED_GAP, MIN_VOTES, SAME_LOCUS, K
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

def semiglobal_dp(read: np.array, window: np.array) -> tuple[np.array, int, int]:
    """Fill the DP matrix (zero first row); return (D, edit_distance, end_col).

    D[i][j] = edit distance between read[:i] and the best substring of window ending at j.
    The zero first row lets the read start anywhere in the window; the minimum of the
    last row lets it end anywhere.
    """

    m, n = len(read), len(window)
    cols = np.arange(n + 1, dtype=np.int32)

    D = np.empty((m + 1, n + 1), dtype=np.int32)
    D[0] = 0

    # Substitution cost for every (read base, window base) pair; N never matches.
    cost = ((read[:, None] != window[None, :]) | (read[:, None] == 4)).astype(np.int32)

    t = np.empty(n + 1, dtype=np.int32)
    for i in range(1, m + 1):
        prev = D[i - 1]
        # Diagonal (match/mismatch) and up (read base vs. gap) terms come from the previous row.
        t[0] = i
        np.minimum(prev[:-1] + cost[i - 1], prev[1:] + 1, out=t[1:])
        # Left term: D[i][j] = min_k (t[k] + j - k), so D[i][j] - j = running min of t[k] - k.
        row = D[i]
        np.subtract(t, cols, out=row)
        np.minimum.accumulate(row, out=row)
        row += cols

    end_col = int(np.argmin(D[m]))
    return D, int(D[m, end_col]), end_col


def traceback(D: np.array, read: np.array, window: np.array, end_col: int) -> int:
    """Return the start column of the alignment ending at end_col.

    Walks back from (len(read), end_col) to row 0, preferring the diagonal move.
    The alignment then covers window[start_col:end_col].
    """

    read, window = read.tolist(), window.tolist()
    i, j = len(read), end_col

    while i > 0:
        d = D[i, j]
        if j > 0:
            r, w = read[i - 1], window[j - 1]
            if d == D[i - 1, j - 1] + (r != w or r == 4):
                i, j = i - 1, j - 1
                continue
        if d == D[i - 1, j] + 1:
            i -= 1
        else:
            j -= 1

    return j


def verify(candidate: tuple[np.array, int], ref_enc: np.array) -> tuple[int, int, int]:
    """Build the reference window, run DP, return (edit_distance, start, end).

    candidate is (read_enc, est_start); start/end are reference coordinates,
    the alignment covers ref_enc[start:end].
    """

    read_enc, est_start = candidate

    lo = max(0, est_start - MARGIN)
    hi = min(len(ref_enc), est_start + len(read_enc) + MARGIN)
    window = ref_enc[lo:hi]

    D, edit_distance, end_col = semiglobal_dp(read_enc, window)
    start_col = traceback(D, read_enc, window, end_col)

    return edit_distance, lo + start_col, lo + end_col


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


def map_read(read_enc: np.array, index: tuple[np.array, np.array], ref_enc: np.array) -> tuple[int, int] | None:
    """Map one read; return (start, end) in the reference or None if unmapped.

    Seeds -> diagonal voting -> DP on the candidates with at least MIN_VOTES hits.
    The best alignment is accepted only if its error rate is at most MAX_ERR_RATE
    and it beats any other location by at least MIN_ED_GAP edits.
    """

    candidates = best_diagonals(*find_hits(read_enc, index))
    candidates = [c for c in candidates if c[1] >= MIN_VOTES]
    if not candidates:
        return None

    # (edit_distance, start, end) for each candidate, best alignment first.
    alignments = sorted(verify((read_enc, est_start), ref_enc) for est_start, _ in candidates)
    edit_distance, start, end = alignments[0]

    if edit_distance > MAX_ERR_RATE * len(read_enc):
        return None

    # A runner-up landing on the same spot comes from hits spilling past the voting window,
    # not from a second copy, so it does not make the read ambiguous.
    if len(alignments) > 1:
        ed2, start2, _ = alignments[1]
        if abs(start2 - start) > SAME_LOCUS and ed2 - edit_distance < MIN_ED_GAP:
            return None

    return start, end

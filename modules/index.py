"""Module for k-mer index construction."""

import numpy as np


def hash_kmers(enc_seq: str, k: int) -> tuple[np.array, np.array]:
    """
    Function to hash the kmers into unique codes and return a boolean mask representing
    valid kmers (not containing N)
    """

    # n is the number of possible k-mers
    n = len(enc_seq) - k + 1

    # Edge cases
    if n <= 0:
        return np.empty(0, dtype=np.uint32), np.empty(0, dtype=bool)

    # Each base gets exactly 2 bits. N (=4) needs 3 bits and would corrupt the neighbouring base so map it to 0.
    bases = (enc_seq % 4).astype(np.uint32)

    # Initialize the codes array for every possible kmer in sequence 
    # Position i in codes represents the k-mer position in sequence.
    codes = np.zeros(n, dtype=np.uint32)
    valid = np.ones(n, dtype=bool)

    # Iterate over all possible kmer positionsi n every kmer
    for j in range(k):
        # Every iteration change the values of the codes based on current slice
        codes = codes * 4 + bases[j:j + n]
        valid = valid & (enc_seq[j:j + n] != 4)

    return codes, valid


def build_index(ref_enc: np.array, k:int) -> tuple[np.array, np.array]:
    """Return (sorted_codes, positions) for all valid reference k-mers."""

    codes, valid = hash_kmers(ref_enc, k)

    positions = np.nonzero(valid)[0].astype(np.uint32)
    codes = codes[valid]

    # Sort by code so identical k-mers are next to each other (stable keeps positions ascending)
    order = np.argsort(codes, kind="stable")
    sorted_codes = codes[order]
    positions = positions[order]

    return sorted_codes, positions

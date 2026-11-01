"""Module for k-mer index construction."""

import numpy as np


def hash_kmers(enc_seq, k):
    """
    Function to hash the kmers into unique codes and return a boolean mask representing
    valid kmers (not containing N)
    """

    # n is the number of possible k-mers
    n = len(enc_seq) - k + 1

    # Edge cases
    if n <= 0:
        return np.empty(0, dtype=np.uint32), np.empty(0, dtype=bool)

    # N (=4) would spill into a third bit so we map it to 0.
    bases = (enc_seq % 4).astype(np.uint32)

    # Append base j of every k-mer at once: shift left by 2 bits, OR in the new base.
    codes = np.zeros(n, dtype=np.uint32)
    valid = np.ones(n, dtype=bool)

    # Iterate over 
    for j in range(k):
        # Every iteration change the values of the codes based on current slice
        codes = codes * 4 + bases[j:j + n]
        valid = valid & (enc_seq[j:j + n] != 4)

    return codes, valid


def build_index(ref_enc, k):
    """Return (sorted_codes, positions) for all valid reference k-mers."""


    raise NotImplementedError
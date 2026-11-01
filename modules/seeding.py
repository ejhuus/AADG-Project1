"""Step 5 — seeding with exact k-mer hits."""

import numpy as np

from modules.config import K, MAX_OCC
from modules.index import kmer_codes


def find_hits(read_enc, index):
    """Return arrays (ref_pos, read_pos) of all exact seed hits.

    k-mers with more than MAX_OCC occurrences in the reference are ignored.
    """
    raise NotImplementedError

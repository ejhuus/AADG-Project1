"""Step 6 — candidate selection by diagonal voting."""

import numpy as np

from modules.config import BAND


def best_diagonals(ref_pos, read_pos):
    """Return the top candidate windows as a list of (est_start, votes).

    The first entry is the densest BAND-wide window of diagonals, the second
    (if any) is the best non-overlapping runner-up.
    """
    raise NotImplementedError

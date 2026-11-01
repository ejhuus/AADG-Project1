"""Step 7 — verification with semi-global edit-distance DP."""

import numpy as np

from modules.config import MARGIN


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

"""Steps 4–8 — map a single read."""

from modules.config import MAX_ERR_RATE, MIN_ED_GAP, MIN_VOTES
from modules.parse_encode import revcomp
from modules.seeding import find_hits
from modules.voting import best_diagonals
from modules.alignment import verify


def map_read(read_enc, index, ref_enc):
    """Map one read (both strands); return (start, end) or None if unmapped."""
    raise NotImplementedError

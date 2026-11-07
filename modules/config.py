"""Tunable parameters (see README, "Parameters")."""

K = 15               # k-mer (seed) length,
MAX_OCC = 30         # ignore k-mers in reference occurring more often than this
BAND = 80            # width of the diagonal voting window (bp)
MIN_VOTES = 10       # minimum seed hits to consider a location
MARGIN = 120         # extra reference on each side of the DP window (bp)
MAX_ERR_RATE = 0.16  # maximum accepted edit distance / read length
MIN_ED_GAP = 10      # required edit-distance gap between best and runner-up

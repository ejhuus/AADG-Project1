"""Module responsible for parsing the FASTA files and encoding the string into numpy array (A=0, C=1, G=2, T=3, other=4)."""

from pathlib import Path
from Bio import SeqIO
import numpy as np

# Lookup table: ASCII byte -> base code. Everything not ACGT (e.g. N) maps to 4.

# Initialize and fill the table with 4 so every other character maps to 4.
LUT = np.full(256, 4, dtype=np.uint8)
for code, bases in enumerate(("Aa", "Cc", "Gg", "Tt")):
    for b in bases:
        LUT[ord(b)] = code


def read_reference(path: Path) -> str:
    """Read the reference sequence file and return a single string."""

    fasta = SeqIO.read(path, 'fasta')
    sequence = str(fasta.seq)

    return sequence


def read_fasta(path: Path):
    """Read sequences to be mapped into a biopython generator yielding read_id and seq."""

    for record in SeqIO.parse(path, "fasta"):
        yield record.id, str(record.seq)


def encode_str(seq: str) -> np.ndarray:
    """Encode as string into a numpy array, using a 256-entry lookup table."""

    raw = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)

    return LUT[raw]


def revcomp(enc):
    """Reverse complement in encoded space (3 - x, reversed; N stays N)."""
    raise NotImplementedError


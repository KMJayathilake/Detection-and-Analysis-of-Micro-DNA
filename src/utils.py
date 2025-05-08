#helper

__all__ = ["revcomp"]

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (A/C/G/T only)."""
    return seq.translate(_COMPLEMENT)[::-1]

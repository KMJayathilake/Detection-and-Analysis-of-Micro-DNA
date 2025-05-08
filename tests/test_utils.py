from src.utils import revcomp

def test_revcomp_basic():
    assert revcomp("ATGC") == "GCAT"
    assert revcomp("acgt") == "acgt"[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))

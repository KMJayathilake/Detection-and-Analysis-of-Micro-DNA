import pytest
import pysam
from pathlib import Path
from src.score import score_circles


CIRCLES = [
    {
        "chrom": "chr1",
        "start": 1400,
        "end": 1500,
        "size_bp": 100,
        "left_tags": {"TAGG"},
        "right_tags": {"TAGG"},
        "support": 3,
    }
]

def mock_mean_depth(*args, **kwargs):
    
    start = args[2]
    return 20 if start == 1400 else 10

def make_empty_bam(path: Path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as outbam:
        pass

def test_score_basic(monkeypatch, tmp_path):
    monkeypatch.setattr("src.score._mean_depth", mock_mean_depth, raising=True)

    dummy_bam = tmp_path / "dummy.bam"
    make_empty_bam(dummy_bam)

    scored = score_circles(
        CIRCLES,
        bam_path=dummy_bam,
        flank=50,
        min_support=1,
        min_score=0.0,
        w_support=1.0,
        w_tag=1.0,
        w_depth=0.5,
    )
    assert len(scored) == 1
    
    assert pytest.approx(scored[0]["score"], rel=1e-2) == 3.585

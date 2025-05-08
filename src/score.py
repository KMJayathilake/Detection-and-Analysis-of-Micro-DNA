#  compute scores for candidate microDNA circles

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List

import pysam

# helper


def _mean_depth(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, max_depth: int = 800) -> float:
    
    total = 0
    length = end - start
    for pileup in bam.pileup(chrom, start, end, truncate=True, max_depth=max_depth):
        total += pileup.nsegments
    return total / length if length else 0.0




def score_circles(
    circles: List[Dict[str, object]],
    bam_path: str | Path,
    *,
    flank: int = 200,
    w_support: float = 1.0,
    w_tag: float = 1.0,
    w_depth: float = 0.5,
    min_support: int = 2,
    min_score: float = 2.0,
) -> List[Dict[str, object]]:
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    scored: List[Dict[str, object]] = []

    for c in circles:
        support = int(c["support"]) if c.get("support") else 1
        
        uniq_left = len(c["left_tags"])
        uniq_right = len(c["right_tags"])
        tag_agreement = 1.0 / max(uniq_left, uniq_right)

        
        chrom = c["chrom"]
        start = int(c["start"])
        end = int(c["end"])
        inside_depth = _mean_depth(bam, chrom, start, end)
        left_bg = _mean_depth(bam, chrom, max(0, start - flank), start)
        right_bg = _mean_depth(bam, chrom, end, end + flank)
        flank_depth = (left_bg + right_bg) / 2 or 1  
        depth_ratio = inside_depth / flank_depth

        
        score = (
            w_support * math.log2(support)
            + w_tag * tag_agreement
            + w_depth * depth_ratio
        )

        c.update(
            {
                "support": support,
                "tag_agreement": round(tag_agreement, 3),
                "depth_ratio": round(depth_ratio, 3),
                "score": round(score, 3),
            }
        )

        if support >= min_support and score >= min_score:
            scored.append(c)

    bam.close()
    return scored



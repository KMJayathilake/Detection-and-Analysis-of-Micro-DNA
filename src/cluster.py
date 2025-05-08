# pair left/right soft clips into candidate microDNA circles

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List

from .utils import revcomp

Circle = Dict[str, object]


def cluster_junctions(
    junction_reads: List[Dict[str, object]],
    *,
    max_gap: int = 1000,
    allow_revcomp: bool = True,
) -> List[Circle]:
    
    # split by chrom & side for quick access
    by_chrom: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: {"left": [], "right": []})
    for rec in junction_reads:
        by_chrom[rec["chrom"]][rec["side"]].append(rec)

    # sort each side records by position for fast scanning
    for chrom in by_chrom:
        by_chrom[chrom]["left"].sort(key=lambda r: r["pos"])
        by_chrom[chrom]["right"].sort(key=lambda r: r["pos"])

    circles: Dict[tuple, Circle] = {}

    # iterate per chromosome
    for chrom, sides in by_chrom.items():
        left_reads = sides["left"]
        right_reads = sides["right"]
        if not left_reads or not right_reads:
            continue

        # two‑pointer sweep
        right_idx = 0
        for L in left_reads:
           
            # move right_idx to first right_read with pos >= L.pos - max_gap

            while right_idx < len(right_reads) and right_reads[right_idx]["pos"] < L["pos"] - max_gap:
                right_idx += 1
            j = right_idx
            while j < len(right_reads) and right_reads[j]["pos"] <= L["pos"]:
                R = right_reads[j]
                # dist check
                if L["pos"] - R["pos"] > max_gap:
                    j += 1
                    continue

                # sequence match
                if L["clip_seq"] == R["clip_seq"] or (
                    allow_revcomp and L["clip_seq"] == revcomp(R["clip_seq"])
                ):
                    start, end = R["pos"], L["pos"]  # start < end
                    key = (chrom, start, end)
                    if key not in circles:
                        circles[key] = {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "size_bp": end - start,
                            "left_tags": set(),
                            "right_tags": set(),
                            "support": 0,
                        }
                    circles[key]["left_tags"].add(L["clip_seq"])
                    circles[key]["right_tags"].add(R["clip_seq"])
                    circles[key]["support"] += 1
                j += 1

    return list(circles.values())

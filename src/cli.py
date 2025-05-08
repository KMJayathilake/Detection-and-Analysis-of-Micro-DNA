#pipline for the project

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import List

import click

from .junctions import find_soft_clipped_reads
from .cluster import cluster_junctions
from .score import score_circles
from .parallel import find_junctions_parallel

logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO, stream=sys.stderr)
logger = logging.getLogger("microdna.cli")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--bam", required=True, type=click.Path(exists=True, dir_okay=False), help="Coordinate sorted BAM/CRAM file")
@click.option("--out", default="circles.tsv", show_default=True, type=click.Path(dir_okay=False), help="Output TSV for final circles")
@click.option("--min-clip", default=10, show_default=True, help="Minimum soft‑clip length")
@click.option("--max-gap", default=1000, show_default=True, help="Max genomic distance between left/right clips (bp)")
@click.option("--flank", default=200, show_default=True, help="Flank size (bp) for depth background")
@click.option("--min-support", default=2, show_default=True, help="Minimum read pair support to keep circle")
@click.option("--min-score", default=2.0, show_default=True, help="Minimum confidence score to report circle")
@click.option("--procs", default=1, show_default=True, help="CPU processes for junction scan (>=1)")
def main(
    bam: str,
    out: str,
    min_clip: int,
    max_gap: int,
    flank: int,
    min_support: int,
    min_score: float,
    procs: int,
):
    #Detect microDNA circles and write TSV report

    bam_path = Path(bam)
    if not bam_path.exists():
        raise SystemExit(f"Input BAM/CRAM not found: {bam_path}")

    
    # 1) soft‑clip detection

    logger.info("Step 1/3  scanning for soft‑clipped reads (%d proc)…", procs)
    if procs > 1:
        junction_reads = find_junctions_parallel(bam_path, min_clip_len=min_clip, procs=procs)
    else:
        junction_reads = find_soft_clipped_reads(bam_path, min_clip_len=min_clip)
    logger.info("  → %d soft‑clipped segments", len(junction_reads))

    
    # 2) clustering

    logger.info("Step 2/3  clustering into candidate circles …")
    circles = cluster_junctions(junction_reads, max_gap=max_gap)
    logger.info("  → %d raw candidate circles", len(circles))

    if not circles:
        logger.warning("No candidate circles found – exiting.")
        return

    
    # 3) scoring + filtering
    
    logger.info("Step 3/3 – scoring & filtering …")
    circles = score_circles(
        circles,
        bam_path,
        flank=flank,
        min_support=min_support,
        min_score=min_score,
    )
    logger.info("  → %d high confidence circles", len(circles))

    if not circles:
        logger.warning("No circles passed filters – nothing written.")
        return

    # write TSV
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = [
        "chrom",
        "start",
        "end",
        "size_bp",
        "support",
        "tag_agreement",
        "depth_ratio",
        "score",
    ]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for c in sorted(circles, key=lambda d: (d["chrom"], d["start"])):
            writer.writerow({k: c.get(k, "-") for k in fieldnames})
    logger.info("Report written → %s", out_path)


if __name__ == "__main__":
    main()

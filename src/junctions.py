# detect soft clipped reads for microDNA pipeline

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pysam

__all__ = ["find_soft_clipped_reads", "JunctionRead"]


logger = logging.getLogger(__name__)

def find_soft_clipped_reads(
    bam_path: str | Path,
    *,
    min_clip_len: int = 10,
    reference_fasta: Optional[str | Path] = None,
    include_secondary: bool = False,
) -> List[JunctionRead]:
    

    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)

    mode = "rc" if bam_path.suffix == ".cram" else "rb"
    bam = pysam.AlignmentFile(
        bam_path,
        mode,
        reference_filename=str(reference_fasta) if reference_fasta else None,
    )

    logger.info("Scanning %s for soft‑clipped reads (min_clip_len=%d)…", bam_path, min_clip_len)

    junction_reads: List[JunctionRead] = []

    for read in bam:
        
        if read.is_unmapped:
            continue
        if not include_secondary and (read.is_secondary or read.is_supplementary):
            continue
        if not read.cigartuples:
            continue  

        
        op_left, len_left = read.cigartuples[0]
        if op_left == 4 and len_left >= min_clip_len:
            junction_reads.append(
                {
                    "chrom": read.reference_name,
                    "pos": read.reference_start,  
                    "clip_seq": read.query_sequence[:len_left],
                    "side": "left",
                    "read_name": read.query_name,
                    "mapq": read.mapping_quality,
                }
            )

        
        op_right, len_right = read.cigartuples[-1]
        if op_right == 4 and len_right >= min_clip_len:
            junction_reads.append(
                {
                    "chrom": read.reference_name,
                    "pos": read.reference_end,  
                    "clip_seq": read.query_sequence[-len_right:],
                    "side": "right",
                    "read_name": read.query_name,
                    "mapq": read.mapping_quality,
                }
            )

    bam.close()
    logger.info("Found %d soft‑clipped segments.", len(junction_reads))
    return junction_reads



# CLI test

if __name__ == "__main__":
    import argparse, csv, sys

    parser = argparse.ArgumentParser(description="Extract soft clipped reads (junctions)")
    parser.add_argument("bam", help="coordinate sorted BAM/CRAM file")
    parser.add_argument("--min", dest="min_clip", type=int, default=10, help="minimum clip length [10]")
    parser.add_argument("--out", default="junctions.tsv", help="output TSV file [junctions.tsv]")
    parser.add_argument("--ref", dest="ref", help="reference fasta (for CRAM)")
    args = parser.parse_args()

    reads = find_soft_clipped_reads(args.bam, min_clip_len=args.min_clip, reference_fasta=args.ref)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["chrom", "pos", "clip_seq", "side", "read_name", "mapq"], delimiter="\t")
        writer.writeheader()
        writer.writerows(reads)

    sys.stderr.write(f"Wrote {len(reads)} rows to {args.out}\n")

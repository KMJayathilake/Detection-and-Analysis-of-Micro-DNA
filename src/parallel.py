#multiprocessing helpers for the microDNA pipeline

from __future__ import annotations

import itertools
import multiprocessing as mp
from pathlib import Path
from typing import List, Dict, Tuple

import pysam

from .junctions import find_soft_clipped_reads  

Junction = Dict[str, object]

__all__ = ["find_junctions_parallel"]



def _contig_list(bam_path: Path) -> List[str]:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return [sq["SN"] for sq in bam.header["SQ"]]


def _process_contig(args: Tuple[Path, str, int]):
    bam_path, contig, min_clip = args
    junctions: List[Junction] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or not read.cigartuples:
                continue
            # left clip
            op, length = read.cigartuples[0]
            if op == 4 and length >= min_clip:
                junctions.append(
                    {
                        "chrom": contig,
                        "pos": read.reference_start,
                        "clip_seq": read.query_sequence[:length],
                        "side": "left",
                        "read_name": read.query_name,
                        "mapq": read.mapping_quality,
                    }
                )
            # right clip
            op, length = read.cigartuples[-1]
            if op == 4 and length >= min_clip:
                junctions.append(
                    {
                        "chrom": contig,
                        "pos": read.reference_end,
                        "clip_seq": read.query_sequence[-length:],
                        "side": "right",
                        "read_name": read.query_name,
                        "mapq": read.mapping_quality,
                    }
                )
    return junctions


# public API


def find_junctions_parallel(
    bam_path: str | Path,
    *,
    min_clip_len: int = 10,
    procs: int = 4,
) -> List[Junction]:
    
    bam_path = Path(bam_path)
    if procs <= 1:
        
        return find_soft_clipped_reads(bam_path, min_clip_len=min_clip_len)

    contigs = _contig_list(bam_path)
    ctx = mp.get_context("spawn")  
    with ctx.Pool(processes=procs) as pool:
        args_iter = ((bam_path, c, min_clip_len) for c in contigs)
        results = pool.map(_process_contig, args_iter)

    
    return list(itertools.chain.from_iterable(results))

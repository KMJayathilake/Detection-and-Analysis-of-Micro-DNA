!!!I had to remove the BAM file due to GitHub’s file size limit. Please add your own BAM file to the data/ folder before running the pipeline.!!!

Python Libraries :
 pysam -	library for reading BAM files, parsing CIGAR strings.
 argparse	-  the CLI script.
 multiprocessing	- Parallel processing across chromosomes.
 csv, os, json, datetime, pathlib	
 pandas.
 matplotlib / samplot

 Tools:
   samtools - BAM indexing and header inspection.
   igv	- Visualization of circle loci and 
   samplot	- Lightweight alternative to IGV for rendering junction read plots.
   conda	


Run the microDNA detection algorithm :
  python -m src.cli \
  --bam data/SRR413984.sorted.NC_000001.10.bam \
  --out output/circles.tsv \
  --min-clip 12 \
  --max-gap 1000 \
  --min-support 2 \
  --min-score 2.0 \
  --procs 8



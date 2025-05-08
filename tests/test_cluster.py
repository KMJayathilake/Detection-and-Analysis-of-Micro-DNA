from src.cluster import cluster_junctions


JUNCTION_READS = [
    {
        "chrom": "chr1",
        "pos": 1500,
        "clip_seq": "TAGG",
        "side": "left",
        "read_name": "readL",
        "mapq": 60,
    },
    {
        "chrom": "chr1",
        "pos": 1400,
        "clip_seq": "TAGG",  
        "side": "right",
        "read_name": "readR",
        "mapq": 60,
    },
]

def test_cluster_single_circle():
    circles = cluster_junctions(JUNCTION_READS, max_gap=200)
    assert len(circles) == 1
    c = circles[0]
    assert (c["chrom"], c["start"], c["end"]) == ("chr1", 1400, 1500)
    assert c["support"] == 1

# output helper for the microdna pipeline

from __future__ import annotations

import csv
import datetime as _dt
import getpass
import json
from pathlib import Path
from typing import Dict, List, Optional

Circle = Dict[str, object]

__all__ = ["write_circle_table", "make_igv_batch"]


# TSV / CSV 


def write_circle_table(
    circles: List[Circle],
    out_path: str | Path,
    *,
    meta: Optional[Dict[str, str]] = None,
    delimiter: str = "\t",
):
   
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "chrom",
        "start",
        "end",
        "size_bp",
        "support",
        "tag_agreement",
        "depth_ratio",
        "score",
    ]
    meta = meta or {}
    meta.update({
        "date": _dt.datetime.now().isoformat(timespec="seconds"),
        "user": getpass.getuser(),
        "pipeline": "microdna-py v0.1",
    })

    with out_path.open("w", newline="") as fh:
        # write metadata as JSON in comment lines (#)
        fh.write("# microDNA circle callset\n")
        for k, v in meta.items():
            fh.write(f"# {k}: {v}\n")
        fh.write("#\n")
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for c in sorted(circles, key=lambda d: (d["chrom"], d["start"])):
            writer.writerow({k: c.get(k, "-") for k in fieldnames})

    return out_path


# IGV helper

def make_igv_batch(
    circles: List[Circle],
    *,
    bam: str | Path,
    genome: str | Path,
    out: str | Path,
    snapshot_dir: Optional[str | Path] = None,
    margin: int = 100,
):
    out = Path(out)
    snapshot_dir = Path(snapshot_dir) if snapshot_dir else out.parent / "snaps"
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    lines = [
        "new",
        f"genome {Path(genome).resolve()}",
        f"load {Path(bam).resolve()}",
        f"snapshotDirectory {snapshot_dir.resolve()}",
    ]

    for idx, c in enumerate(circles, 1):
        chrom, start, end = c["chrom"], int(c["start"]), int(c["end"])
        view_start = max(0, start - margin)
        view_end = end + margin
        lines.append(f"goto {chrom}:{view_start}-{view_end}")
        snap_name = f"circle_{idx}_{chrom}_{start}_{end}.png"
        lines.append(f"snapshot {snap_name}")

    lines.append("exit")

    out.write_text("\n".join(lines) + "\n")
    return out



def write_html_summary(circles: List[Circle], out_html: str | Path):
    """Very small HTML table; good for quick sharing."""
    out_html = Path(out_html)
    rows = "\n".join(
        f"<tr><td>{c['chrom']}:{c['start']}-{c['end']}</td><td>{c['size_bp']}</td><td>{c['support']}</td><td>{c['score']:.2f}</td></tr>"
        for c in circles
    )
    html = f"""
    <html><head><title>microDNA circles</title></head><body>
    <h2>High-confidence microDNA circles</h2>
    <table border="1" cellspacing="0" cellpadding="4">
      <tr><th>Location</th><th>Size (bp)</th><th>Support</th><th>Score</th></tr>
      {rows}
    </table>
    </body></html>"""
    out_html.write_text(html)
    return out_html

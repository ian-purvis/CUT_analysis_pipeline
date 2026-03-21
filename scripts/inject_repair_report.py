#!/usr/bin/env python3
"""
Inject paired-read repair stats into FastQC HTML report.
Adds a "Paired Read Repair" section showing reads before/after and orphaned reads removed.
"""
import argparse
import gzip
import re
import sys
from pathlib import Path


def count_reads_gzip(path: Path) -> int:
    """Count reads in gzipped FASTQ (lines / 4)."""
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            return sum(1 for _ in f) // 4
    except (OSError, ValueError):
        return -1


def parse_repair_log(log_path: Path) -> dict:
    """Parse repair.sh log for stats. Returns dict with keys like read1, read2, pairs, orphaned_r1, orphaned_r2."""
    stats = {}
    if not log_path.exists():
        return stats
    text = log_path.read_text(errors="ignore")
    # Common bbmap output patterns
    for pattern, key in [
        (r"Read\s+1:\s*(\d+)", "read1"),
        (r"Read\s+2:\s*(\d+)", "read2"),
        (r"Pairs:\s*(\d+)", "pairs"),
        (r"Orphaned\s+R1:\s*(\d+)", "orphaned_r1"),
        (r"Orphaned\s+R2:\s*(\d+)", "orphaned_r2"),
    ]:
        m = re.search(pattern, text, re.I)
        if m:
            stats[key] = int(m.group(1))
    return stats


def build_repair_html(
    sample: str,
    before_r1: int,
    before_r2: int,
    after_r1: int,
    after_r2: int,
    log_stats: dict,
) -> str:
    """Build HTML section for repair report."""
    pairs_kept = min(after_r1, after_r2) if after_r1 >= 0 and after_r2 >= 0 else 0
    orphan_r1 = log_stats.get("orphaned_r1")
    orphan_r2 = log_stats.get("orphaned_r2")
    if orphan_r1 is None and before_r1 >= 0 and pairs_kept >= 0:
        orphan_r1 = max(0, before_r1 - pairs_kept)
    if orphan_r2 is None and before_r2 >= 0 and pairs_kept >= 0:
        orphan_r2 = max(0, before_r2 - pairs_kept)

    total_removed = (orphan_r1 or 0) + (orphan_r2 or 0)
    changed = total_removed > 0

    rows = []
    rows.append(("<tr><td>Reads before repair (R1)</td><td>{:,}</td></tr>").format(before_r1 if before_r1 >= 0 else 0))
    rows.append(("<tr><td>Reads before repair (R2)</td><td>{:,}</td></tr>").format(before_r2 if before_r2 >= 0 else 0))
    rows.append(("<tr><td>Paired reads kept</td><td>{:,}</td></tr>").format(pairs_kept))
    if orphan_r1 is not None and orphan_r1 > 0:
        rows.append(("<tr><td>Orphaned R1 removed</td><td>{:,}</td></tr>").format(orphan_r1))
    if orphan_r2 is not None and orphan_r2 > 0:
        rows.append(("<tr><td>Orphaned R2 removed</td><td>{:,}</td></tr>").format(orphan_r2))
    if total_removed > 0:
        rows.append(("<tr><td><strong>Total reads removed</strong></td><td><strong>{:,}</strong></td></tr>").format(total_removed))

    status_text = "Reads were repaired (orphaned reads removed)" if changed else "No repair needed (reads already paired)"
    status_style = "color:#F0AD4E" if changed else "color:#5CB85C"  # warning orange / success green

    html = f"""
<div class="module">
<h2>Paired Read Repair</h2>
<p>bbmap repair.sh was run to fix paired FASTQs with unequal read counts. Trim Galore requires equal R1/R2 counts.</p>
<table>
<thead><tr><th>Metric</th><th>Count</th></tr></thead>
<tbody>
{"".join(rows)}
</tbody>
</table>
<p><strong>Status:</strong> <span style="{status_style}">{status_text}</span></p>
<p><em>Orphaned reads</em> have no matching mate in the other file (e.g. truncated run, dropped reads). They are excluded so downstream paired-end tools receive properly synchronized files.</p>
</div>
"""
    return html


def inject_into_fastqc(html_path: Path, repair_html: str) -> None:
    """Inject repair section before first FastQC module."""
    content = html_path.read_text(encoding="utf-8", errors="replace")
    marker = '<div class="module">'
    if marker in content:
        content = content.replace(marker, repair_html + marker, 1)
    else:
        # Fallback: insert after <div class="main">
        marker = '<div class="main">'
        if marker in content:
            content = content.replace(marker, marker + repair_html, 1)
        else:
            sys.stderr.write(f"Could not find insertion point in {html_path}\n")
            return
    html_path.write_text(content, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description="Inject repair stats into FastQC HTML")
    ap.add_argument("--sample", required=True, help="Sample name")
    ap.add_argument("--input-r1", required=True, help="Input R1 FASTQ (before repair)")
    ap.add_argument("--input-r2", required=True, help="Input R2 FASTQ (before repair)")
    ap.add_argument("--output-r1", required=True, help="Output R1 FASTQ (after repair)")
    ap.add_argument("--output-r2", required=True, help="Output R2 FASTQ (after repair)")
    ap.add_argument("--repair-log", required=True, help="repair.sh log file")
    ap.add_argument("--fastqc-html", required=True, help="FastQC HTML to modify (R1 report)")
    ap.add_argument("--fastqc-html-r2", default="", help="FastQC HTML for R2 (optional, inject into both)")
    args = ap.parse_args()

    input_r1 = Path(args.input_r1)
    input_r2 = Path(args.input_r2)
    output_r1 = Path(args.output_r1)
    output_r2 = Path(args.output_r2)
    repair_log = Path(args.repair_log)
    fastqc_html = Path(args.fastqc_html)

    before_r1 = count_reads_gzip(input_r1)
    before_r2 = count_reads_gzip(input_r2)
    after_r1 = count_reads_gzip(output_r1)
    after_r2 = count_reads_gzip(output_r2)
    log_stats = parse_repair_log(repair_log)

    repair_html = build_repair_html(
        args.sample, before_r1, before_r2, after_r1, after_r2, log_stats
    )

    if fastqc_html.exists():
        inject_into_fastqc(fastqc_html, repair_html)
    if args.fastqc_html_r2:
        html_r2 = Path(args.fastqc_html_r2)
        if html_r2.exists():
            inject_into_fastqc(html_r2, repair_html)


if __name__ == "__main__":
    main()

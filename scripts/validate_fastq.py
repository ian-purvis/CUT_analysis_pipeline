#!/usr/bin/env python3
"""
Validate FASTQ files for format integrity. Detects corrupted/incomplete files
via gzip integrity check (CRC) and format parsing. Writes valid sample names to output file.
Always exits 0 so the pipeline continues; invalid samples are skipped.
"""
import gzip
import subprocess
import sys
from pathlib import Path


def gzip_integrity_ok(path: Path) -> bool:
    """Run gzip -t to catch CRC errors and truncation. Fast and catches corruption anywhere in file."""
    if path.suffix != ".gz":
        return True  # Not gzipped, skip this check
    try:
        subprocess.run(
            ["gzip", "-t", str(path)],
            capture_output=True,
            check=True,
            timeout=600,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def validate_fastq(path: str, max_records: int = 1800000) -> bool:
    """Check FASTQ format by parsing first N records (~7M lines). Catches truncation and mid-file corruption."""
    path = Path(path)
    if not path.exists():
        return False
    # First: gzip integrity (catches CRC errors, truncation)
    if path.suffix == ".gz" and not gzip_integrity_ok(path):
        return False
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt") as f:
            for i, line in enumerate(f):
                if i >= max_records * 4:
                    break
                line_num = i % 4
                if line_num == 0:
                    if not line.startswith("@"):
                        return False
                elif line_num == 2:
                    if not line.startswith("+"):
                        return False
        return True
    except Exception:
        return False


def main():
    # Usage: validate_fastq.py <fastq_dir> <sample1> <sample2> ... --output <out.txt>
    if len(sys.argv) < 4 or "--output" not in sys.argv:
        print("Usage: validate_fastq.py <fastq_dir> <sample1> <sample2> ... --output <out.txt>", file=sys.stderr)
        sys.exit(2)

    idx = sys.argv.index("--output")
    fastq_dir = Path(sys.argv[1].strip().rstrip("/")).resolve()
    samples = sys.argv[2:idx]
    out_path = sys.argv[idx + 1] if idx + 1 < len(sys.argv) else "samples.txt"

    valid_samples = []
    for i, sample in enumerate(samples):
        print(f"[validate_fastq] Checking {sample} ({i+1}/{len(samples)})...", file=sys.stderr)
        r1 = fastq_dir / sample / f"{sample}_R1.fastq.gz"
        r2 = fastq_dir / sample / f"{sample}_R2.fastq.gz"
        if validate_fastq(str(r1)) and validate_fastq(str(r2)):
            valid_samples.append(sample)
        else:
            print(f"[validate_fastq] Skipping {sample}: corrupted or incomplete FASTQ", file=sys.stderr)

    out = Path(out_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\n".join(valid_samples) + "\n")

    print(f"[validate_fastq] Valid samples ({len(valid_samples)}/{len(samples)}): {', '.join(valid_samples)}", file=sys.stderr)
    sys.exit(0)  # Always succeed - we've written the valid list


if __name__ == "__main__":
    main()

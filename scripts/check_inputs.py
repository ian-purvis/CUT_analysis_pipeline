#!/usr/bin/env python3
"""
Input validation for CUT&RUN pipeline modules.
Each check_* function validates inputs for a module and exits with 1 on failure.
Usage: check_inputs.py <module> [args...]
"""
import os
import sys
from pathlib import Path


def get_samples_from_file(samples_file: str) -> list:
    """Read validated sample list from file."""
    p = Path(samples_file)
    if not p.exists():
        return []
    with open(p) as f:
        return [s.strip() for s in f if s.strip()]


def get_samples_from_config(config_path: str) -> list:
    """Read sample list from config (before validation). Fallback when discovery finds nothing."""
    try:
        import yaml
    except ImportError:
        return []
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    samples = []
    for reps in cfg.get("samples", {}).values():
        samples.extend(reps)
    return samples


def discover_samples(fastq_dir: Path) -> list:
    """Discover samples from FASTQ directory: subdirs with {sample}_R1.fastq.gz and {sample}_R2.fastq.gz."""
    samples = []
    if not fastq_dir.is_dir():
        return samples
    for sub in sorted(fastq_dir.iterdir()):
        if not sub.is_dir():
            continue
        name = sub.name
        r1 = sub / f"{name}_R1.fastq.gz"
        r2 = sub / f"{name}_R2.fastq.gz"
        if r1.exists() and r2.exists():
            samples.append(name)
    return samples


def check_qc(fastq_dir: str, config_path: str, workdir: str, samples_out: str = None) -> list:
    """Validate inputs for QC module: discover samples from FASTQ dir, ensure R1/R2 exist."""
    fastq_dir = fastq_dir.strip().rstrip("/")  # Normalize to avoid path drift
    fastq_path = (Path(workdir) / fastq_dir).resolve()
    if not fastq_path.exists():
        print(f"ERROR: FASTQ directory not found: {fastq_path}", file=sys.stderr)
        sys.exit(1)
    samples = discover_samples(fastq_path)
    if not samples:
        samples = get_samples_from_config(config_path)
    if not samples:
        print("ERROR: No samples found. Expected subdirs with {sample}_R1.fastq.gz and {sample}_R2.fastq.gz", file=sys.stderr)
        sys.exit(1)
    missing = []
    for s in samples:
        r1 = fastq_path / s / f"{s}_R1.fastq.gz"
        r2 = fastq_path / s / f"{s}_R2.fastq.gz"
        if not r1.exists():
            missing.append(f"{s}_R1.fastq.gz")
        if not r2.exists():
            missing.append(f"{s}_R2.fastq.gz")
    if missing:
        print(f"ERROR: Missing FASTQ files: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    if samples_out:
        out_path = Path(samples_out).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            f.write("\n".join(samples) + "\n")
    print(f"OK: {len(samples)} samples: {', '.join(samples)}", file=sys.stderr)
    return samples


def check_align(genome_bwa: str, samples_file: str, trimmed_dir: str, workdir: str) -> list:
    """Validate inputs for alignment: genome, BWA index, trimmed FASTQs."""
    genome = Path(genome_bwa)
    if not genome.is_absolute():
        genome = Path(workdir) / genome_bwa
    if not genome.exists():
        print(f"ERROR: Genome not found: {genome}", file=sys.stderr)
        sys.exit(1)
    # BWA index: .bwt or .ann
    idx = genome.with_suffix(genome.suffix + ".bwt")
    if not idx.exists():
        idx = genome.with_suffix(genome.suffix + ".ann")
    if not idx.exists():
        print(f"ERROR: BWA index not found. Run: bwa index {genome}", file=sys.stderr)
        sys.exit(1)
    samples = get_samples_from_file(samples_file)
    if not samples:
        print(f"ERROR: No validated samples. Run QC module first to create {samples_file}", file=sys.stderr)
        sys.exit(1)
    trimmed = Path(workdir) / trimmed_dir
    missing = []
    for s in samples:
        r1 = trimmed / f"{s}_R1_val_1.fq.gz"
        if not r1.exists():
            missing.append(r1.name)
    if missing:
        print(f"ERROR: Missing trimmed FASTQs. Run QC module first. Missing: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    print(f"OK: Genome, BWA index, and {len(samples)} trimmed samples ready", file=sys.stderr)
    return samples


def check_dedup(samples_file: str, bam_dir: str, workdir: str) -> list:
    """Validate inputs for deduplication: BAM exists for each sample."""
    samples = get_samples_from_file(samples_file)
    if not samples:
        print(f"ERROR: No validated samples. Run QC module first.", file=sys.stderr)
        sys.exit(1)
    bam_path = Path(workdir) / bam_dir
    missing = []
    for s in samples:
        bam = bam_path / f"{s}.bam"
        if not bam.exists():
            missing.append(bam.name)
    if missing:
        print(f"ERROR: Missing BAM files. Run alignment module first. Missing: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    print(f"OK: {len(samples)} BAM files ready for deduplication", file=sys.stderr)
    return samples


def check_peaks(samples_file: str, bam_dir: str, control: str, workdir: str) -> list:
    """Validate inputs for peak calling: control + treatment dedup BAMs."""
    samples = get_samples_from_file(samples_file)
    if not samples:
        print(f"ERROR: No validated samples.", file=sys.stderr)
        sys.exit(1)
    tf_samples = [s for s in samples if s != control]
    if not tf_samples:
        print(f"ERROR: No TF samples (all are control?)", file=sys.stderr)
        sys.exit(1)
    bam_path = Path(workdir) / bam_dir
    ctrl_bam = bam_path / f"{control}.dedup.bam"
    if not ctrl_bam.exists():
        print(f"ERROR: Control BAM not found: {ctrl_bam}. Run dedup module first.", file=sys.stderr)
        sys.exit(1)
    missing = []
    for s in tf_samples:
        bam = bam_path / f"{s}.dedup.bam"
        if not bam.exists():
            missing.append(bam.name)
    if missing:
        print(f"ERROR: Missing dedup BAMs for peak calling. Run dedup module first. Missing: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    print(f"OK: Control {control} and {len(tf_samples)} TF samples ready for peak calling", file=sys.stderr)
    return tf_samples


def check_bigwig(samples_file: str, bam_dir: str, workdir: str) -> list:
    """Validate inputs for BigWig: dedup BAM files exist."""
    samples = get_samples_from_file(samples_file)
    if not samples:
        print(f"ERROR: No validated samples.", file=sys.stderr)
        sys.exit(1)
    bam_path = Path(workdir) / bam_dir
    missing = []
    for s in samples:
        bam = bam_path / f"{s}.dedup.bam"
        if not bam.exists():
            missing.append(bam.name)
    if missing:
        print(f"ERROR: Missing dedup BAMs. Run dedup module first. Missing: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    print(f"OK: {len(samples)} dedup BAM files ready for BigWig", file=sys.stderr)
    return samples


def main():
    if len(sys.argv) < 2:
        print("Usage: check_inputs.py <module> [args...]", file=sys.stderr)
        sys.exit(2)
    module = sys.argv[1].lower()
    workdir = os.getcwd()
    if len(sys.argv) > 2 and sys.argv[2] == "--workdir":
        workdir = sys.argv[3]
        args = sys.argv[4:]
    else:
        args = sys.argv[2:]

    if module == "qc":
        if len(args) < 2:
            print("Usage: check_inputs.py qc <fastq_dir> <config.yaml> [samples_out]", file=sys.stderr)
            sys.exit(2)
        samples_out = args[2] if len(args) > 2 else None
        check_qc(args[0], args[1], workdir, samples_out)
    elif module == "align":
        if len(args) < 3:
            print("Usage: check_inputs.py align <genome_bwa> <samples.txt> <trimmed_dir>", file=sys.stderr)
            sys.exit(2)
        check_align(args[0], args[1], args[2], workdir)
    elif module == "dedup":
        if len(args) < 2:
            print("Usage: check_inputs.py dedup <samples.txt> <bam_dir>", file=sys.stderr)
            sys.exit(2)
        check_dedup(args[0], args[1], workdir)
    elif module == "peaks":
        if len(args) < 3:
            print("Usage: check_inputs.py peaks <samples.txt> <bam_dir> <control>", file=sys.stderr)
            sys.exit(2)
        check_peaks(args[0], args[1], args[2], workdir)
    elif module == "bigwig":
        if len(args) < 2:
            print("Usage: check_inputs.py bigwig <samples.txt> <bam_dir>", file=sys.stderr)
            sys.exit(2)
        check_bigwig(args[0], args[1], workdir)
    else:
        print(f"Unknown module: {module}. Use: qc, align, dedup, peaks, bigwig", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

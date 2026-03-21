#!/usr/bin/env bash
# Run alignment module only: BWA mem
# Requires: QC module complete (trimmed FASTQs, samples.txt)
# Usage: ./run_align.sh <fastq_dir> [jobs]
# Resumable: re-run after power/interrupt; skips completed jobs, reruns incomplete.

if [[ -z "$1" ]]; then
  echo "Usage: $0 <fastq_dir> [jobs]"
  exit 1
fi
FASTQ_DIR="${1%/}"   # Remove trailing slash for path consistency
JOBS=${2:-4}
OPTS="-j $JOBS -k --latency-wait 60 --rerun-incomplete"
echo "=== Module 2: Alignment (BWA) ==="
echo "  fastq_dir: $FASTQ_DIR"
snakemake --unlock 2>/dev/null || true
snakemake all_align --config fastq_dir="$FASTQ_DIR" $OPTS

#!/usr/bin/env bash
# Run deduplication module only: samtools markdup
# Requires: Alignment module complete (BAM files)
# Usage: ./run_dedup.sh <fastq_dir> [jobs]
# Resumable: re-run after power/interrupt; skips completed jobs, reruns incomplete.

if [[ -z "$1" ]]; then
  echo "Usage: $0 <fastq_dir> [jobs]"
  exit 1
fi
FASTQ_DIR="${1%/}"   # Remove trailing slash for path consistency
JOBS=${2:-4}
OPTS="-j $JOBS -k --latency-wait 60 --rerun-incomplete"
echo "=== Module 3: Deduplication ==="
echo "  fastq_dir: $FASTQ_DIR"
snakemake --unlock 2>/dev/null || true
snakemake all_dedup --config fastq_dir="$FASTQ_DIR" $OPTS

#!/usr/bin/env bash
# Run peak calling module only: MACS3
# Requires: Dedup module complete (dedup BAMs for control + treatments)
# Usage: ./run_peaks.sh <fastq_dir> [jobs]
# Resumable: re-run after power/interrupt; skips completed jobs, reruns incomplete.

if [[ -z "$1" ]]; then
  echo "Usage: $0 <fastq_dir> [jobs]"
  exit 1
fi
FASTQ_DIR="${1%/}"   # Remove trailing slash for path consistency
JOBS=${2:-4}
OPTS="-j $JOBS -k --latency-wait 60 --rerun-incomplete"
echo "=== Module 4: Peak calling (MACS3) ==="
echo "  fastq_dir: $FASTQ_DIR"
snakemake --unlock 2>/dev/null || true
snakemake all_peaks --config fastq_dir="$FASTQ_DIR" $OPTS

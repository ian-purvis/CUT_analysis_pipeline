#!/usr/bin/env bash
# Run QC module only: setup, validate, FastQC, Trim Galore
# Usage: ./run_qc.sh <fastq_dir> [jobs]
#   fastq_dir: path to FASTQ folder (relative to cutrun_pipeline, or absolute)
#   jobs: parallel jobs (default 4)
# Resumable: re-run after power/interrupt; skips completed jobs, reruns incomplete.

if [[ -z "$1" ]]; then
  echo "Usage: $0 <fastq_dir> [jobs]"
  echo "  fastq_dir: path to folder with sample subfolders (e.g. ../../Mu CUT&Tag_OTX2_2023/FASTQ)"
  echo "  jobs: parallel jobs (default 4)"
  exit 1
fi
FASTQ_DIR="${1%/}"   # Remove trailing slash for path consistency
JOBS=${2:-4}
OPTS="-j $JOBS -k --latency-wait 60 --rerun-incomplete"
echo "=== Module 1: QC (setup, validate, FastQC, Trim Galore) ==="
echo "  fastq_dir: $FASTQ_DIR"
snakemake --unlock 2>/dev/null || true
snakemake all_qc --config fastq_dir="$FASTQ_DIR" $OPTS

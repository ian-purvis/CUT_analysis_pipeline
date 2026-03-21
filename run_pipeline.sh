#!/usr/bin/env bash
# Run full CUT&RUN pipeline (all modules in sequence)
# Usage: ./run_pipeline.sh <fastq_dir> [jobs]
#   fastq_dir: path to FASTQ folder (relative to cutrun_pipeline, or absolute)
#   jobs: parallel jobs (default 4)
#
# Runs each module one at a time: QC -> Align -> Dedup -> Peaks -> BigWig
# Each module validates its inputs before running.
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

snakemake --unlock 2>/dev/null || true
echo "=== Full pipeline (all modules) ==="
echo "  fastq_dir: $FASTQ_DIR"
for target in all_qc all_align all_dedup all_peaks all_bigwig; do
    echo ""
    echo ">>> Running $target"
    snakemake $target --config fastq_dir="$FASTQ_DIR" $OPTS || exit 1
done
echo ""
echo "=== Pipeline complete ==="

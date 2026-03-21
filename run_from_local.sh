#!/usr/bin/env bash
# Run pipeline with Snakemake metadata on local WSL filesystem (avoids OneDrive I/O issues).
# Data (FASTQ, results) stays on OneDrive; only .snakemake/ lives in ~/cutrun_work.
#
# Usage: ./run_from_local.sh <fastq_dir> [jobs] [target]
#   fastq_dir: path to FASTQ folder (relative to cutrun_pipeline, or absolute)
#   jobs: parallel jobs (default 4)
#   target: all_qc | all_align | all_dedup | all_peaks | all_bigwig | all (default: all)
#
# Example: ./run_from_local.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4 all

if [[ -z "$1" ]]; then
  echo "Usage: $0 <fastq_dir> [jobs] [target]"
  echo "  fastq_dir: path to folder with sample subfolders"
  echo "  jobs: parallel jobs (default 4)"
  echo "  target: all_qc | all_align | all_dedup | all_peaks | all_bigwig | all (default: all)"
  echo ""
  echo "Runs Snakemake from ~/cutrun_work so .snakemake metadata is local."
  echo "FASTQ and results remain on OneDrive."
  exit 1
fi

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
FASTQ_DIR="${1%/}"
# Resolve to absolute path so Snakemake (run from ~/cutrun_work) uses correct OneDrive locations
if [[ "$FASTQ_DIR" != /* ]]; then
  FASTQ_DIR="$(cd "$PIPELINE_DIR" && cd "$FASTQ_DIR" && pwd)"
fi
RESULTS_DIR="$(dirname "$FASTQ_DIR")/results"
JOBS=${2:-4}
TARGET=${3:-all}
OPTS="-j $JOBS -k --latency-wait 60 --rerun-incomplete"
META_DIR="$HOME/cutrun_work"

mkdir -p "$META_DIR"
cd "$META_DIR" || exit 1

echo "=== Run from local (OneDrive-safe) ==="
echo "  Metadata dir: $META_DIR"
echo "  Pipeline:     $PIPELINE_DIR"
echo "  fastq_dir:    $FASTQ_DIR"
echo "  results_dir:  $RESULTS_DIR"
echo "  target:       $TARGET"
echo ""

snakemake -s "$PIPELINE_DIR/Snakefile" --configfile "$PIPELINE_DIR/config.yaml" \
  --config fastq_dir="$FASTQ_DIR" results_dir="$RESULTS_DIR" --unlock 2>/dev/null || true

if [[ "$TARGET" == "all" ]]; then
  for t in all_qc all_align all_dedup all_peaks all_bigwig; do
    echo ">>> Running $t"
    snakemake -s "$PIPELINE_DIR/Snakefile" --configfile "$PIPELINE_DIR/config.yaml" \
      --config fastq_dir="$FASTQ_DIR" results_dir="$RESULTS_DIR" $OPTS "$t" || exit 1
  done
else
  snakemake -s "$PIPELINE_DIR/Snakefile" --configfile "$PIPELINE_DIR/config.yaml" \
    --config fastq_dir="$FASTQ_DIR" results_dir="$RESULTS_DIR" $OPTS "$TARGET"
fi

echo ""
echo "=== Done ==="

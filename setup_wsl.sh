#!/bin/bash
# Run this script in WSL (Ubuntu) to set up the CUT&RUN pipeline
# The pipeline requires Linux - bioconda packages (bwa, samtools, macs3) don't support Windows natively.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v conda &>/dev/null; then
    echo "Conda not found. Install Miniconda first:"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh"
    echo ""
    echo "Then restart your terminal and run this script again."
    exit 1
fi

echo "=== Creating conda environment (cutrun) ==="
conda env create -f environment.yaml

echo ""
echo "=== Setup complete ==="
echo "Activate with: conda activate cutrun"
echo "Then run: snakemake -j 4"
echo ""
echo "To index mm39 genome first: ./index_mm39.sh"

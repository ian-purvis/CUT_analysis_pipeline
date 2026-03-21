#!/bin/bash
# Setup mm39 (GRCm39) reference genome for CUT&RUN pipeline
# Run from cutrun_pipeline/ directory: bash setup_reference.sh

set -e
REF_DIR="reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "Downloading mm39 reference genome from UCSC..."
if [ ! -f mm39.fa ]; then
    if [ ! -f mm39.fa.gz ]; then
        wget -q --show-progress "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" || \
        curl -L -o mm39.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
    fi
    echo "Decompressing..."
    gunzip -k mm39.fa.gz 2>/dev/null || gunzip mm39.fa.gz
fi

echo "Building BWA index (this may take 10-20 minutes)..."
bwa index mm39.fa

echo "Building samtools index..."
samtools faidx mm39.fa

echo "mm39 reference setup complete: $(pwd)/mm39.fa"

# Optional: mm39 blacklist (excluderanges)
# Download from excluderanges if available
echo ""
echo "Optional: mm39 blacklist for BigWig filtering"
echo "You can create reference/mm39-blacklist.bed manually from excluderanges R package,"
echo "or leave blacklist_bed empty in config.yaml to skip blacklist filtering."

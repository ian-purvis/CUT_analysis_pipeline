#!/bin/bash
# Build BWA and samtools indexes for mm39 genome (run in WSL with conda env activated)
# Reads genome path from config.yaml - edit config.yaml for your system

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Extract genome_bwa from config.yaml (handles quoted paths)
GENOME=$(grep 'genome_bwa:' config.yaml | sed 's/.*genome_bwa:[[:space:]]*"\([^"]*\)".*/\1/')
if [[ -z "$GENOME" ]]; then
    echo "ERROR: Could not read genome_bwa from config.yaml"
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: Genome not found at $GENOME"
    echo "Edit config.yaml to set the correct genome path, then copy the genome file."
    exit 1
fi

echo "Building samtools index..."
samtools faidx "$GENOME"

echo "Building BWA index (~10-20 min)..."
bwa index "$GENOME"

echo "Done."

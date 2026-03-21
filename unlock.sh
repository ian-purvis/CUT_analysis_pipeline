#!/usr/bin/env bash
# Clear Snakemake lock after power loss or forced kill.
# Run this if you get "Error: Directory cannot be locked" or similar.
# Usage: ./unlock.sh

snakemake --unlock

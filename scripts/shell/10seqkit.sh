#!/bin/bash
#SBATCH -J filter_contigs
#SBATCH -o filter_contigs.%j.out
#SBATCH -e filter_contigs.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00

# ============================================================
# Filter metagenomic contigs by length using SeqKit
#
# Samples: a1-a120
#
# Input:
#   metaSPAdes scaffold assemblies
#
# Output:
#   Contigs >= 2,000 bp for downstream Plasmer analysis
#
# Before running:
#   1. Make sure SeqKit is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
# ============================================================

set -euo pipefail

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/assembly/metaSPAdes"
OUTPUT_DIR="${PROJECT_DIR}/metagenome/assembly/filtered_contigs"

MIN_LENGTH=2000

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    INPUT_FILE="${INPUT_DIR}/${SAMPLE}_scaffolds.fasta"
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}_filtered.fasta"

    echo ">>> Filtering contigs for ${SAMPLE}..."

    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "WARNING: Missing input file: ${INPUT_FILE}"
        continue
    fi

    seqkit seq \
        -m "${MIN_LENGTH}" \
        -g \
        "${INPUT_FILE}" \
        -o "${OUTPUT_FILE}"

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> Contig filtering completed for all samples."
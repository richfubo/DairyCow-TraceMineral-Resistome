#!/bin/bash
#SBATCH -J Plasmer
#SBATCH -o Plasmer.%j.out
#SBATCH -e Plasmer.%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=480G
#SBATCH --time=72:00:00

# ============================================================
# Plasmid-like and chromosome-like contig classification
# using Plasmer
#
# Shotgun metagenomic data from dairy cow fecal samples
#
# Samples: a1-a120
#
# Input:
#   Filtered metagenomic contigs generated from metaSPAdes
#
# Output:
#   Plasmer predictions used for genomic-context analysis
#   of ARGs, MRGs, and MGEs
#
# Before running:
#   1. Make sure Plasmer is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
#   3. Set PLASMER_DB to the local Plasmer database.
# ============================================================

set -euo pipefail

THREADS=16

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/assembly/filtered_contigs"
OUTPUT_DIR="${PROJECT_DIR}/metagenome/Plasmer"

PLASMER_DB="/path/to/Plasmer/database"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    INPUT_FASTA="${INPUT_DIR}/${SAMPLE}_filtered.fasta"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"

    echo ">>> Running Plasmer for ${SAMPLE}..."

    if [[ ! -f "${INPUT_FASTA}" ]]; then
        echo "WARNING: Missing input file: ${INPUT_FASTA}"
        continue
    fi

    mkdir -p "${SAMPLE_OUT}"

    Plasmer \
        -g "${INPUT_FASTA}" \
        -p "${SAMPLE}" \
        -d "${PLASMER_DB}" \
        -t "${THREADS}" \
        -o "${SAMPLE_OUT}"

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> Plasmer analysis completed for all samples."
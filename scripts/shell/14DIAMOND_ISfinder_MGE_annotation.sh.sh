#!/bin/bash
#SBATCH -J ISfinder_annotation
#SBATCH -o ISfinder_annotation.%j.out
#SBATCH -e ISfinder_annotation.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=300:00:00

# ============================================================
# Mobile genetic element annotation against the ISfinder
# database using DIAMOND
#
# Samples: a1-a120
#
# Input:
#   Protein sequences predicted by Prodigal from
#   metagenomic contigs
#
# Output:
#   DIAMOND blastp results against the ISfinder database
#
# Before running:
#   1. Make sure DIAMOND is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
#   3. Set ISFINDER_DB to the local DIAMOND-formatted
#      ISfinder database.
# ============================================================

set -euo pipefail

# --------------------------
# Parameters
# --------------------------
THREADS=8

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

# Prodigal output
INPUT_DIR="${PROJECT_DIR}/metagenome/Prodigal"

# ISfinder annotation output
OUTPUT_DIR="${PROJECT_DIR}/metagenome/ISfinder"

# DIAMOND-formatted ISfinder database
ISFINDER_DB="/path/to/ISfinder/ISfinder"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    FAA="${INPUT_DIR}/${SAMPLE}/${SAMPLE}.proteins.faa"
    OUT="${OUTPUT_DIR}/${SAMPLE}.isfinder.tsv"

    if [[ ! -f "${FAA}" ]]; then
        echo "WARNING: Missing protein file: ${FAA}"
        continue
    fi

    echo ">>> Running ISfinder annotation for ${SAMPLE}..."

    diamond blastp \
        -q "${FAA}" \
        -d "${ISFINDER_DB}" \
        -o "${OUT}" \
        -f 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --id 60 \
        --query-cover 70 \
        --evalue 1e-5 \
        --max-target-seqs 10 \
        --threads "${THREADS}"

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> ISfinder annotation completed for all samples."
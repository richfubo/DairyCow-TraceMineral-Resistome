#!/bin/bash
#SBATCH -J BacMet_annotation
#SBATCH -o BacMet_annotation.%j.out
#SBATCH -e BacMet_annotation.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=300:00:00

# ============================================================
# Metal resistance gene annotation against the BacMet2
# database using DIAMOND
#
# Samples: a1-a120
#
# Input:
#   Protein sequences predicted by Prodigal from
#   metagenomic contigs
#
# Output:
#   DIAMOND blastp results against the BacMet2 database
#
# Before running:
#   1. Make sure DIAMOND is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
#   3. Set BACMET_DB to the local DIAMOND-formatted
#      BacMet2 database.
# ============================================================

set -euo pipefail

# --------------------------
# Parameters
# --------------------------
THREADS=8

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

# Prodigal output
INPUT_DIR="${PROJECT_DIR}/metagenome/Prodigal"

# BacMet2 annotation output
OUTPUT_DIR="${PROJECT_DIR}/metagenome/BacMet2"

# DIAMOND-formatted BacMet2 database
BACMET_DB="/path/to/BacMet2/BacMet2"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    FAA="${INPUT_DIR}/${SAMPLE}/${SAMPLE}.proteins.faa"
    OUT="${OUTPUT_DIR}/${SAMPLE}.bacmet.tsv"

    if [[ ! -f "${FAA}" ]]; then
        echo "WARNING: Missing protein file: ${FAA}"
        continue
    fi

    echo ">>> Running BacMet2 annotation for ${SAMPLE}..."

    diamond blastp \
        -q "${FAA}" \
        -d "${BACMET_DB}" \
        -o "${OUT}" \
        -f 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --id 70 \
        --query-cover 70 \
        --evalue 1e-10 \
        --threads "${THREADS}"

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> BacMet2 annotation completed for all samples."
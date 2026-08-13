#!/bin/bash
#SBATCH -J SARG_annotation
#SBATCH -o SARG_annotation.%j.out
#SBATCH -e SARG_annotation.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=300:00:00

# ============================================================
# ARG annotation against the SARG database using DIAMOND
#
# Samples: a1-a120
#
# Input:
#   Protein sequences predicted by Prodigal from
#   metagenomic contigs
#
# Output:
#   DIAMOND blastp results against the SARG database
#
# Before running:
#   1. Make sure DIAMOND is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
#   3. Set SARG_DB to the local DIAMOND-formatted SARG database.
# ============================================================

set -euo pipefail

# --------------------------
# Parameters
# --------------------------
THREADS=8

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

# Prodigal output
INPUT_DIR="${PROJECT_DIR}/metagenome/Prodigal"

# SARG annotation output
OUTPUT_DIR="${PROJECT_DIR}/metagenome/SARG"

# DIAMOND-formatted SARG database
SARG_DB="/path/to/SARG/sarg"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    FAA="${INPUT_DIR}/${SAMPLE}/${SAMPLE}.proteins.faa"
    OUT="${OUTPUT_DIR}/${SAMPLE}.sarg.tsv"

    if [[ ! -f "${FAA}" ]]; then
        echo "WARNING: Missing protein file: ${FAA}"
        continue
    fi

    echo ">>> Running SARG annotation for ${SAMPLE}..."

    diamond blastp \
        -q "${FAA}" \
        -d "${SARG_DB}" \
        -o "${OUT}" \
        -f 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --id 70 \
        --query-cover 70 \
        --evalue 1e-10 \
        --threads "${THREADS}"

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> SARG annotation completed for all samples."
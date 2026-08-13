#!/bin/bash
#SBATCH -J BacMet_stage2
#SBATCH -o BacMet_stage2.%j.out
#SBATCH -e BacMet_stage2.%j.err
#SBATCH --cpus-per-task=20
#SBATCH --time=300:00:00

# ============================================================
# BacMet Stage Two analysis
# Shotgun metagenomic data from dairy cow fecal samples
#
# Samples: a1-a120
#
# Input:
#   BacMet Stage One output
#
# Required files:
#   extracted.fa
#   meta_data_online.txt
#
# Output:
#   BacMet Stage Two results for downstream MRG profiling
#
# Before running:
#   1. Modify PROJECT_DIR according to the local environment
#   2. Set BACMET_PIPELINE to the local installation of the
#      BacMet Stage Two pipeline
# ============================================================

set -euo pipefail

# --------------------------
# Parameters
# --------------------------
THREADS=20

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/BacMet/stage_one"
OUTPUT_DIR="${PROJECT_DIR}/metagenome/BacMet/stage_two"

# Custom BacMet Stage Two pipeline
BACMET_PIPELINE="/path/to/argoap_pipeline_stagetwo_version2_BacMet"

# --------------------------
# Input files
# --------------------------
EXTRACTED_FASTA="${INPUT_DIR}/extracted.fa"
METADATA="${INPUT_DIR}/meta_data_online.txt"
METADATA_NUMERIC="${INPUT_DIR}/meta_data_online.numeric.txt"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Check required input files
# --------------------------
if [[ ! -f "${EXTRACTED_FASTA}" ]]; then
    echo "ERROR: ${EXTRACTED_FASTA} not found."
    exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: ${METADATA} not found."
    exit 1
fi

# --------------------------
# Convert the first metadata column
# to sequential numeric IDs
# --------------------------
echo ">>> Preparing numeric metadata..."

awk 'BEGIN{FS=OFS="\t"}
     NR==1 {print; next}
     {$1=NR-1; print}' \
     "${METADATA}" > "${METADATA_NUMERIC}"

# --------------------------
# Run BacMet Stage Two
# --------------------------
echo ">>> Starting BacMet Stage Two..."
echo "Samples         : a1-a120"
echo "Input directory : ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"

"${BACMET_PIPELINE}" \
    -i "${EXTRACTED_FASTA}" \
    -m "${METADATA_NUMERIC}" \
    -o "${OUTPUT_DIR}" \
    -n "${THREADS}"

echo ">>> BacMet Stage Two completed successfully."
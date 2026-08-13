#!/bin/bash
#SBATCH -J kneaddata_batch
#SBATCH -o kneaddata_batch.%j.out
#SBATCH -e kneaddata_batch.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=72:00:00

# ============================================================
# Host read removal and quality filtering using KneadData
# Shotgun metagenomic sequencing data from dairy cow feces
#
# Samples: a1-a120
#
# Input:
#   a1_clean_R1.fq.gz
#   a1_clean_R2.fq.gz
#   ...
#   a120_clean_R1.fq.gz
#   a120_clean_R2.fq.gz
#
# Host reference:
#   Bos taurus reference genome
#
# Before running:
#   1. Make sure KneadData is available in $PATH
#   2. Modify PROJECT_DIR and DB_PATH for your environment
# ============================================================

# --------------------------
# Parameters
# --------------------------
THREADS=8

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/clean_reads"
OUTPUT_DIR="${PROJECT_DIR}/metagenome/02_kneaddata"

# KneadData database built from the Bos taurus reference genome
DB_PATH="/path/to/Bos_taurus_kneaddata_database"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    R1="${INPUT_DIR}/${SAMPLE}_clean_R1.fq.gz"
    R2="${INPUT_DIR}/${SAMPLE}_clean_R2.fq.gz"

    echo ">>> Running KneadData for ${SAMPLE}..."

    # Check that both paired-end files exist
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "Warning: input files for ${SAMPLE} were not found. Skipping."
        continue
    fi

    kneaddata \
        --input1 "${R1}" \
        --input2 "${R2}" \
        --output "${OUTPUT_DIR}" \
        --output-prefix "${SAMPLE}_kneaddata" \
        --reference-db "${DB_PATH}" \
        --threads "${THREADS}" \
        --log "${OUTPUT_DIR}/${SAMPLE}_kneaddata.log" \
        --remove-intermediate-output \
        --decontaminate-pairs lenient \
        --reorder \
        --bypass-trf

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> KneadData processing completed for all samples."
#!/bin/bash
#SBATCH -J fastqc_multiqc
#SBATCH -o fastqc.%j.out
#SBATCH -e fastqc.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# ============================================================
# Quality assessment of shotgun metagenomic sequencing reads
#
# Samples: a1-a120
# Tools: FastQC and MultiQC
#
# Before running:
#   1. Make sure FastQC and MultiQC are available in $PATH
#   2. Modify PROJECT_DIR according to your local environment
# ============================================================

set -euo pipefail

# --------------------------
# User-defined directories
# --------------------------
PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/clean_reads"
OUT_DIR="${PROJECT_DIR}/metagenome/01_fastqc"

mkdir -p "${OUT_DIR}"

# --------------------------
# Parallel settings
# --------------------------
MAX_JOBS=4
THREADS_PER_JOB=2

echo "Starting FastQC analysis..."

job_count=0

# Samples a1-a120
for ID in {1..120}; do

    SAMPLE="a${ID}"

    for READ in R1 R2; do

        FILE="${INPUT_DIR}/${SAMPLE}_clean_${READ}.fq.gz"

        if [[ ! -f "${FILE}" ]]; then
            echo "Warning: ${FILE} not found. Skipping."
            continue
        fi

        echo "Processing ${FILE} ..."

        fastqc \
            --nogroup \
            -t "${THREADS_PER_JOB}" \
            -o "${OUT_DIR}" \
            "${FILE}" &

        ((job_count+=1))

        if (( job_count % MAX_JOBS == 0 )); then
            wait
        fi

    done
done

wait

echo "Generating MultiQC report..."

multiqc "${OUT_DIR}" \
    -o "${OUT_DIR}"

echo "FastQC and MultiQC analyses completed successfully."
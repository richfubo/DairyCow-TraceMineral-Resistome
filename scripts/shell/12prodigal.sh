#!/bin/bash
#SBATCH -J Prodigal
#SBATCH -o Prodigal.%j.out
#SBATCH -e Prodigal.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --time=300:00:00

# ============================================================
# ORF prediction using Prodigal
# Shotgun metagenomic data from dairy cow fecal samples
#
# Samples: a1-a120
#
# Input:
#   Metagenomic contigs >= 2,000 bp
#
# Output:
#   Predicted protein sequences, nucleotide coding sequences,
#   and GFF annotations for downstream ARG, MRG, and MGE
#   annotation.
#
# Before running:
#   1. Make sure Prodigal is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
# ============================================================

set -euo pipefail

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

INPUT_DIR="${PROJECT_DIR}/metagenome/assembly/filtered_contigs"
OUTPUT_DIR="${PROJECT_DIR}/metagenome/Prodigal"

mkdir -p "${OUTPUT_DIR}"

# --------------------------
# Process samples a1-a120
# --------------------------
for ID in {1..120}; do

    SAMPLE="a${ID}"

    CONTIG="${INPUT_DIR}/${SAMPLE}_filtered.fasta"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"

    if [[ ! -f "${CONTIG}" ]]; then
        echo "WARNING: Missing contig file: ${CONTIG}"
        continue
    fi

    mkdir -p "${SAMPLE_OUT}"

    echo ">>> Running Prodigal for ${SAMPLE}..."

    prodigal \
        -i "${CONTIG}" \
        -a "${SAMPLE_OUT}/${SAMPLE}.proteins.faa" \
        -d "${SAMPLE_OUT}/${SAMPLE}.genes.fna" \
        -o "${SAMPLE_OUT}/${SAMPLE}.genes.gff" \
        -f gff \
        -p meta \
        -q

    echo ">>> Finished ${SAMPLE}"

done

echo ">>> Prodigal ORF prediction completed for all samples."
#!/bin/bash
#SBATCH -J HUMAnN_postprocess
#SBATCH -o HUMAnN_postprocess.%j.out
#SBATCH -e HUMAnN_postprocess.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00

# ============================================================
# Merge and normalize HUMAnN functional profiles
#
# Samples: a1-a120
#
# Input:
#   Per-sample HUMAnN outputs
#
# Output:
#   Merged and normalized MetaCyc pathway abundance tables
#   for downstream functional analyses
#
# Main downstream file:
#   pathabundance_relab_unstratified.tsv
#
# Before running:
#   1. Make sure HUMAnN utilities are available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
# ============================================================

set -euo pipefail

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome"

# Per-sample HUMAnN results
RESULT_DIR="${PROJECT_DIR}/metagenome/HUMAnN"

# Temporary directory used to collect per-sample tables
FLAT_DIR="${PROJECT_DIR}/metagenome/HUMAnN_postprocessing/flat"

# Merged HUMAnN tables
MERGED_DIR="${PROJECT_DIR}/metagenome/HUMAnN_postprocessing/merged"

mkdir -p "${FLAT_DIR}"
mkdir -p "${MERGED_DIR}"

echo ">>> Checking HUMAnN output files for a1-a120..."

missing_count=0

for ID in {1..120}; do

    SAMPLE="a${ID}"
    SAMPLE_DIR="${RESULT_DIR}/${SAMPLE}"

    GF="${SAMPLE_DIR}/${SAMPLE}_merged_genefamilies.tsv"
    PA="${SAMPLE_DIR}/${SAMPLE}_merged_pathabundance.tsv"
    PC="${SAMPLE_DIR}/${SAMPLE}_merged_pathcoverage.tsv"

    if [[ ! -f "${GF}" || ! -f "${PA}" || ! -f "${PC}" ]]; then
        echo "WARNING: Missing HUMAnN output for ${SAMPLE}"
        missing_count=$((missing_count + 1))
    fi

done

if [[ "${missing_count}" -ne 0 ]]; then
    echo "ERROR: ${missing_count} samples have missing HUMAnN outputs."
    exit 1
fi

echo ">>> All 120 samples passed the file check."

# --------------------------
# Collect per-sample HUMAnN tables
# --------------------------
rm -rf "${FLAT_DIR}"
mkdir -p "${FLAT_DIR}"

for ID in {1..120}; do

    SAMPLE="a${ID}"

    ln -sf \
        "${RESULT_DIR}/${SAMPLE}/${SAMPLE}_merged_genefamilies.tsv" \
        "${FLAT_DIR}/"

    ln -sf \
        "${RESULT_DIR}/${SAMPLE}/${SAMPLE}_merged_pathabundance.tsv" \
        "${FLAT_DIR}/"

    ln -sf \
        "${RESULT_DIR}/${SAMPLE}/${SAMPLE}_merged_pathcoverage.tsv" \
        "${FLAT_DIR}/"

done

# --------------------------
# Merge HUMAnN tables
# --------------------------
echo ">>> Merging gene family tables..."

humann_join_tables \
    --input "${FLAT_DIR}" \
    --output "${MERGED_DIR}/genefamilies.tsv" \
    --file_name genefamilies

echo ">>> Merging pathway abundance tables..."

humann_join_tables \
    --input "${FLAT_DIR}" \
    --output "${MERGED_DIR}/pathabundance.tsv" \
    --file_name pathabundance

echo ">>> Merging pathway coverage tables..."

humann_join_tables \
    --input "${FLAT_DIR}" \
    --output "${MERGED_DIR}/pathcoverage.tsv" \
    --file_name pathcoverage

# --------------------------
# Normalize pathway abundance
# to relative abundance
# --------------------------
echo ">>> Renormalizing pathway abundance..."

humann_renorm_table \
    --input "${MERGED_DIR}/pathabundance.tsv" \
    --output "${MERGED_DIR}/pathabundance_relab.tsv" \
    --units relab

# --------------------------
# Split stratified and
# unstratified pathway profiles
# --------------------------
echo ">>> Splitting stratified and unstratified tables..."

humann_split_stratified_table \
    --input "${MERGED_DIR}/pathabundance_relab.tsv" \
    --output "${MERGED_DIR}"

# --------------------------
# Check final output
# --------------------------
UNSTRATIFIED="${MERGED_DIR}/pathabundance_relab_unstratified.tsv"

if [[ ! -f "${UNSTRATIFIED}" ]]; then
    echo "ERROR: Unstratified pathway table was not generated."
    exit 1
fi

COLUMN_NUMBER=$(head -n 1 "${UNSTRATIFIED}" | tr '\t' '\n' | wc -l)

echo ">>> Columns in final table: ${COLUMN_NUMBER}"
echo ">>> Expected: 121 columns (1 pathway + 120 samples)"

if [[ "${COLUMN_NUMBER}" -eq 121 ]]; then
    echo ">>> Sample number check passed."
else
    echo "WARNING: Unexpected number of columns."
fi

echo ">>> HUMAnN postprocessing completed successfully."
echo ">>> Final pathway table:"
echo "${UNSTRATIFIED}"
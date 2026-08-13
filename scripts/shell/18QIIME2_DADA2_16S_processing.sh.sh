#!/bin/bash
#SBATCH -J QIIME2_16S
#SBATCH -o QIIME2_16S.%j.out
#SBATCH -e QIIME2_16S.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# ============================================================
# 16S rRNA gene sequence processing using QIIME 2 and DADA2
#
# Samples: a1-a120
#
# Amplicon region:
#   V3-V4
#
# Workflow:
#   FASTQ import
#      -> DADA2 denoising
#      -> ASV table
#      -> taxonomic assignment
#      -> removal of organelle-associated sequences
#      -> export for downstream statistical analyses in R
#
# Before running:
#   1. Make sure QIIME 2 is available in $PATH.
#   2. Modify PROJECT_DIR according to the local environment.
#   3. Set CLASSIFIER to the appropriate SILVA 138 classifier.
# ============================================================

set -euo pipefail

THREADS=12

PROJECT_DIR="/path/to/DairyCow-TraceMineral-Resistome/16S"

MANIFEST="${PROJECT_DIR}/metadata/manifest.tsv"
METADATA="${PROJECT_DIR}/metadata/metadata.tsv"

WORK_DIR="${PROJECT_DIR}/qiime2"
EXPORT_DIR="${PROJECT_DIR}/exported"

# SILVA 138 classifier
CLASSIFIER="/path/to/SILVA138_classifier.qza"

mkdir -p "${WORK_DIR}"
mkdir -p "${EXPORT_DIR}"

cd "${WORK_DIR}"

# ------------------------------------------------------------
# Step 1. Import paired-end reads
# ------------------------------------------------------------
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${MANIFEST}" \
    --output-path demux-paired-end.qza \
    --input-format PairedEndFastqManifestPhred33V2

# ------------------------------------------------------------
# Step 2. Summarize sequencing quality
# ------------------------------------------------------------
qiime demux summarize \
    --i-data demux-paired-end.qza \
    --o-visualization demux-summary.qzv

# ------------------------------------------------------------
# Step 3. DADA2 denoising
# ------------------------------------------------------------
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux-paired-end.qza \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 222 \
    --p-trunc-len-r 229 \
    --p-n-threads "${THREADS}" \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza

# ------------------------------------------------------------
# Step 4. Summarize DADA2 output
# ------------------------------------------------------------
qiime metadata tabulate \
    --m-input-file denoising-stats.qza \
    --o-visualization denoising-stats.qzv

qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file "${METADATA}"

qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv

# ------------------------------------------------------------
# Step 5. Taxonomic assignment using SILVA 138
# ------------------------------------------------------------
qiime feature-classifier classify-sklearn \
    --i-classifier "${CLASSIFIER}" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

# ------------------------------------------------------------
# Step 6. Remove mitochondrial and chloroplast sequences
# ------------------------------------------------------------
qiime taxa filter-table \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table table-no-organelle.qza

qiime taxa filter-seqs \
    --i-sequences rep-seqs.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-sequences rep-seqs-no-organelle.qza

# ------------------------------------------------------------
# Step 7. Remove singleton ASVs
# ------------------------------------------------------------
qiime feature-table filter-features \
    --i-table table-no-organelle.qza \
    --p-min-frequency 2 \
    --o-filtered-table table-final.qza

# ------------------------------------------------------------
# Step 8. Generate taxonomic composition visualization
# ------------------------------------------------------------
qiime taxa barplot \
    --i-table table-final.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file "${METADATA}" \
    --o-visualization taxa-bar-plots.qzv

# ------------------------------------------------------------
# Step 9. Export ASV table, taxonomy, and representative
# sequences for downstream analyses in R
# ------------------------------------------------------------
qiime tools export \
    --input-path table-final.qza \
    --output-path "${EXPORT_DIR}/feature-table"

qiime tools export \
    --input-path taxonomy.qza \
    --output-path "${EXPORT_DIR}/taxonomy"

qiime tools export \
    --input-path rep-seqs-no-organelle.qza \
    --output-path "${EXPORT_DIR}/representative-sequences"

echo ">>> QIIME 2 processing completed successfully."
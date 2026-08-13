#!/bin/bash

#SBATCH -J fastqc_limited_parallel
#SBATCH -o /public/home/2021180/job/metaanalysis/fastqc.%j.out
#SBATCH -e /public/home/2021180/job/metaanalysis/fastqc.%j.err
#SBATCH --partition=Cnode_all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8      


export PATH=/public/home/2021180/miniconda3/envs/fastqc/bin:$PATH
export PATH=/public/home/2021180/miniconda3/envs/multiqc/bin:$PATH

INPUT_DIR="/public/home/2021180/users/fubo/data/gut_microbiome_cow/megegenomics_MbPL202308519/2_cleandata"
OUT_DIR="/public/home/2021180/users/fubo/data/gut_microbiome_cow/metaresult/01fastqc"
mkdir -p "$OUT_DIR"


MAX_JOBS=4        
THREADS_PER_JOB=2 

FILES=($(ls "$INPUT_DIR"/*.fq.gz))

echo "Start FastQC with limited concurrency..."

job_count=0
for FILE in "${FILES[@]}"; do
    echo "Processing $FILE ..."
    fastqc --nogroup -t $THREADS_PER_JOB -o "$OUT_DIR" "$FILE" &

    ((job_count++))

   
    if (( job_count % MAX_JOBS == 0 )); then
        wait
    fi
done


wait


echo "Generating MultiQC report..."
multiqc "$OUT_DIR" -o "$OUT_DIR"

echo "All done ✅"

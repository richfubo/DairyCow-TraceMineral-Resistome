#!/bin/bash

#SBATCH -J fastqc_limited_parallel
#SBATCH -o /public/home/2021180/job/metaanalysis/fastqc.%j.out
#SBATCH -e /public/home/2021180/job/metaanalysis/fastqc.%j.err
#SBATCH --partition=Cnode_all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8      # 总共分配 8 核，根据并发数合理设定

# 加载环境
export PATH=/public/home/2021180/miniconda3/envs/fastqc/bin:$PATH
export PATH=/public/home/2021180/miniconda3/envs/multiqc/bin:$PATH

INPUT_DIR="/public/home/2021180/users/fubo/data/gut_microbiome_cow/megegenomics_MbPL202308519/2_cleandata"
OUT_DIR="/public/home/2021180/users/fubo/data/gut_microbiome_cow/metaresult/01fastqc"
mkdir -p "$OUT_DIR"

# 控制最大并发数
MAX_JOBS=4        # 同时最多跑 4 个 fastqc
THREADS_PER_JOB=2 # 每个 fastqc 使用 2 个线程

FILES=($(ls "$INPUT_DIR"/*.fq.gz))

echo "Start FastQC with limited concurrency..."

job_count=0
for FILE in "${FILES[@]}"; do
    echo "Processing $FILE ..."
    fastqc --nogroup -t $THREADS_PER_JOB -o "$OUT_DIR" "$FILE" &

    ((job_count++))

    # 如果达到最大并发数，等待全部结束再继续
    if (( job_count % MAX_JOBS == 0 )); then
        wait
    fi
done

# 等待所有剩余任务
wait

# MultiQC 汇总
echo "Generating MultiQC report..."
multiqc "$OUT_DIR" -o "$OUT_DIR"

echo "All done ✅"

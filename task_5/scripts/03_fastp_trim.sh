#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=../logs/fastp_%j.out
#SBATCH --error=../logs/fastp_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00

DATA_DIR=../data
TRIM_DIR=../results/trimmed

mkdir -p $TRIM_DIR

for fq in $DATA_DIR/*.fastq.gz; do
    base=$(basename $fq .fastq.gz)
    fastp \
        -i $fq \
        -o $TRIM_DIR/${base}_trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --cut_right \
        --cut_right_window_size 5 \
        --cut_right_mean_quality 20 \
        --length_required 36 \
        -j $TRIM_DIR/${base}_fastp.json \
        -h $TRIM_DIR/${base}_fastp.html \
        -w 4
done

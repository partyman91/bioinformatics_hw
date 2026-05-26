#!/bin/bash
#SBATCH --job-name=fastqc_after
#SBATCH --output=../logs/fastqc_after_%j.out
#SBATCH --error=../logs/fastqc_after_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00

TRIM_DIR=../results/trimmed
OUT_DIR=../results/fastqc_after
MULTIQC_DIR=../results/multiqc_after

mkdir -p $OUT_DIR $MULTIQC_DIR

fastqc -t 4 -o $OUT_DIR $TRIM_DIR/*_trimmed.fastq.gz

multiqc $OUT_DIR $TRIM_DIR -o $MULTIQC_DIR -n multiqc_after

#!/bin/bash
#SBATCH --job-name=fastqc_before
#SBATCH --output=../logs/fastqc_before_%j.out
#SBATCH --error=../logs/fastqc_before_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00

DATA_DIR=../data
OUT_DIR=../results/fastqc_before
MULTIQC_DIR=../results/multiqc_before

mkdir -p $OUT_DIR $MULTIQC_DIR

fastqc -t 4 -o $OUT_DIR $DATA_DIR/*.fastq.gz

multiqc $OUT_DIR -o $MULTIQC_DIR -n multiqc_before

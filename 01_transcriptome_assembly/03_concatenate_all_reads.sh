#!/bin/bash
#SBATCH -n 20
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH -N 1
#SBATCH --mem=100G


# All the trimmed RNA-Seq data was concatenated seperately for each species for Trinity assembly input.
# Scripts was ran from each species data directory seperately.

cat *_R1.fastq > reads_all_R1.fq
cat *_R2.fastq > reads_all_R2.fq

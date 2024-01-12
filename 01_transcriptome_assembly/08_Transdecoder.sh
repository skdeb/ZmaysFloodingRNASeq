#!/bin/bash
#SBATCH -n 10
#SBATCH -p main
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

module load bio/transdecoder/5.5.0

TransDecoder.LongOrfs -t Trinity_filtered.fasta

TransDecoder.Predict -t Trinity_filtered.fasta --single_best_only

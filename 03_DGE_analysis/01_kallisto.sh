#!/bin/bash
#SBATCH -n 16
#SBATCH -p threaded
#SBATCH -q threaded
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

source ~/.bash_profile
conda activate rnaseq_kallisto
export PATH=$PATH:/mrm/bin/parallel-20211022/:/mrm/bin/parallel-20211022/bin:/mrm/bin/parallel-20211022/share

#The script was run seperately for each species in their respective directories.  

cdna_fa=$1
index_name=$2
reads_dir=$3

kallisto index -i $index_name $cdna_fa &>index.log

# ids is a txt file contains the unique part of the name of each samples.

cat ids | parallel -j 1 "kallisto quant -i $index_name \
 -o ./{} -t 16 --rf-stranded \
 $reads_dir/{}_R1.fq $reads_dir/{}_R2.fq &> {}.log"

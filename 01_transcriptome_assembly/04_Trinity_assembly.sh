#!/bin/bash
#SBATCH -n 10
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu


#The script was run seperately for each species in their respective directories. 

module load singularity/3.7.2

singularity exec -B ${PWD}:${PWD} \
 -e /grps2/mrmckain/bin/trinity_singularity/trinityrnaseq.v2.15.1.simg Trinity \
 --seqType fq \
 --left data/reads_all_R1.fq \
 --right data/reads_all_R2.fq \
 --SS_lib_type RF --max_memory 50G --CPU 10 \
 --output trinity_out &> trinity_out.log



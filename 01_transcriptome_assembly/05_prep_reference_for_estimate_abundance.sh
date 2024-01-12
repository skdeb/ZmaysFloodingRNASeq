#!/bin/bash
#SBATCH -n 20
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

module load singularity/3.7.2
module load bio/bowtie/2.3

singularity exec -B ${PWD}:${PWD} \
-e /grps2/mrmckain/bin/trinity_singularity/trinityrnaseq.v2.15.1.simg \
 /usr/local/bin/util/align_and_estimate_abundance.pl \
 --thread_count 20 --transcripts Trinity.fasta \
 --est_method RSEM --aln_method bowtie2 \
 --trinity_mode --prep_reference

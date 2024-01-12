#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16 
#SBATCH --mem=100G
#SBATCH -p threaded
#SBATCH -q threaded
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu


module load singularity/3.7.2
module load bio/bowtie/2.3
export PATH=$PATH:/mrm/bin/parallel-20211022/:/mrm/bin/parallel-20211022/bin:/mrm/bin/parallel-20211022/share
module load perl/5.22.1threaded

cat ids | parallel -j 1 "singularity exec -B ${PWD}:${PWD} \
 -B /grps2/mrmckain/bin/rnaseq_tools/RSEM-1.3.3:/RSEM \
 -B /scratch/skdeb/flooding_data_analysis/data/trimmed_data/:/data \
 -e /grps2/mrmckain/bin/trinity_singularity/trinityrnaseq.v2.15.1.simg \
 /usr/local/bin/util/align_and_estimate_abundance.pl \
 --thread_count 20 --transcripts Trinity.fasta \
 --seqType fq \
 --left /data/{}_R1.fastq \
 --right /data/{}_R2.fastq \
 --SS_lib_type RF --est_method RSEM --aln_method bowtie2 \
 --trinity_mode --output_dir {} &> {}.log"

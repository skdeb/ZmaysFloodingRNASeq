#!/bin/bash
#SBATCH -n 32
#SBATCH -p ultrahigh
#SBATCH --qos mrmckain
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

export DK_ROOT=/share/apps/dotkit
. /share/apps/dotkit/bash/.dk_init
me=`whoami`

module load bio/fastqc/0.11.5

reads_dir=/scratch/skdeb/flooding_data_analysis/data
#make directory for this step

mkdir 01_fastqc
cd 01_fastqc

#make output directories for each species

mkdir trip_fastqc vossia_fastqc zea_fastqc znic_fastqc

fastqc $reads_dir/trip_data/*R?.fastq -t 32 -o ./trip_fastqc
fastqc $reads_dir/vossia_data/*R?.fastq -t 32 -o ./vossia_fastqc
fastqc $reads_dir/zea_data/*R?.fastq -t 32 -o ./zea_fastqc
fastqc $reads_dir/znic_data/*R?.fastq -t 32 -o ./znic_fastqc

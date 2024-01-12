#!/bin/bash
#SBATCH -n 20
#SBATCH -p main
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

module load singularity/3.7.2

#############################################################
# database only need to be created once. if you need to update the database run step 1.
# this commands will create the database (myTrinotate.sqlite/ or whatever name you give for the database) in your current working directory. You are going to use this database in downstream analysis.
#You will also find a boilerplate.sqlite database within the TRINOTATE_DATA_DIR. The TRINOTATE_DATA_DIR and boilerplate.sqlite can be reused for future Trinotate runs, in which case you can just copy and rename the boilerplate.sqlite database in a new working directory instead of having to rerun --create for future Trinotate runs that would leverage the same set of database resources.

# I have the TRINOTATE_DATA_DIR in here: /grps2/mrmckain/bin/annotation_tools/Functional_annotation/trinotate 
############################################################
# 1. Sequence Databases
singularity exec -B ${PWD}:/anno \
-e /grps2/mrmckain/bin/annotation_tools/Functional_annotation/trinotate/trinotate.v4.0.2.simg \
 /usr/local/src/Trinotate/Trinotate \
 --create \
 --db myTrinotate.sqlite \
 --trinotate_data_dir /anno/TRINOTATE_DATA_DIR \
 --use_diamond


##############################################################

# 2. set path for required files

container_path=/grps2/mrmckain/bin/annotation_tools/Functional_annotation/trinotate/trinotate.v4.0.2.simg
scripts_path=/usr/local/src/Trinotate
db=$1 # keep in the working directory. I copy and rename the boilerplate.sqlite to myTrinotate.sqlite
transcripts=$2  #your target transcriptome in fasta format
peps=$3 #coding regions translated in fasta format
maps=$4 #pairwise mappings between gene and transcript isoform identifiers(created in trinity run)

################################################################

# 3. Inititalize your Trinotate sqlite database with your sequence data

singularity exec -B ${PWD}:/anno \
-e $container_path \
 $scripts_path/Trinotate\
 --db /anno/$db --init \
 --gene_trans_map /anno/$maps \
 --transcript_fasta /anno/$transcripts \
 --transdecoder_pep /anno/$peps

##################################################################

# 4. Running Sequence Analyses

singularity exec -B ${PWD}:/anno \
 -B /grps2/mrmckain/bin/annotation_tools/Functional_annotation/trinotate:/data \
 -e $container_path \
 $scripts_path/Trinotate\
 --db /anno/$db \
 --CPU 10 \
 --transcript_fasta /anno/$transcripts \
 --transdecoder_pep /anno/$peps \
 --trinotate_data_dir /data/TRINOTATE_DATA_DIR \
 --run "swissprot_blastp swissprot_blastx pfam infernal EggnogMapper" \
 --use_diamond \

####################################################################

# 5. Generating the Trinotate Report

singularity exec -B ${PWD}:/anno \
 -e $container_path \
 $scripts_path/Trinotate\
 --db /anno/$db \
 --report > myTrinotate.xls

#####################################################################

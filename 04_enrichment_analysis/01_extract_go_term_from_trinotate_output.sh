#!/bin/bash
#SBATCH -n 10
#SBATCH -p threaded
#SBATCH -q threaded
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH -o slurm_output.%J
#SBATCH -e slurm_error.%J
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skdeb@uahpc.ua.edu

module load singularity/3.7.2
dir=$(pwd)
singularity exec -B ${PWD}:/anno \
-e /grps2/mrmckain/bin/annotation_tools/Functional_annotation/trinotate/trinotate.v4.0.2.simg \
 /usr/local/src/Trinotate/util/extract_GO_assignments_from_Trinotate_xls.pl \
 --Trinotate_xls /anno/Trinotate.xls \
 -G --include_ancestral_terms > go_annotations.txt

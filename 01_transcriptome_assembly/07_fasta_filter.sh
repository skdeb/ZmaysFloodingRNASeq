#!/bin/bash
#SBATCH -n 10
#SBATCH -p main
#SBATCH -N 1
#SBATCH --mem=20G


#filter_fasta_by_rsem_values.pl script is from Trinity utility scripts available on Github.

perl filter_fasta_by_rsem_values.pl \
 --rsem_output=Leaf0_1/RSEM.isoforms.results,Leaf0_2/RSEM.isoforms.results,Leaf0_3/RSEM.isoforms.results,Leaf4_1/RSEM.isoforms.results,Leaf4_2/RSEM.isoforms.results,Leaf4_3/RSEM.isoforms.results,Leaf72_1/RSEM.isoforms.results,Leaf72_2/RSEM.isoforms.results,Leaf72_3/RSEM.isoforms.results,Root0_1/RSEM.isoforms.results,Root0_2/RSEM.isoforms.results,Root0_3/RSEM.isoforms.results,Root4_1/RSEM.isoforms.results,Root4_2/RSEM.isoforms.results,Root4_3/RSEM.isoforms.results,Root72_1/RSEM.isoforms.results,Root72_2/RSEM.isoforms.results,Root72_3/RSEM.isoforms.results \
 --fasta=Trinity.fasta \
 --output=Trinity_filtered.fasta \
 --tpm_cutoff=2 \
 --isopct_cutoff=20

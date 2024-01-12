#!/bin/bash

## GOATOOLs for GO enrichment analysis of diffentially expressed genes
## This scirpt was used for different set of genes at various stages of the analysis.

source ~/.bash_profile
conda activate goatools

study=$1 # list of differentially expressed genes (aka study genes)
population=$2 # list of all present transcripts in the study (aka population)
go_asso=$3 # txt file with genes and their associate go terms (see example file format)
outfile=$4 # output file name txt or any other format
obo=$5
find_enrichment.py --pval=0.05 --method=fdr_bh \
 --pval_field=fdr_bh --obo=$obo\
 $study $population $go_asso \
 --outfile=$outfile

import sys
import pandas as pd

'''
Editing the GO annotations format to use with GOATOOLS.
GO-term annotation data of Zea mays was downloaded in csv format from
here: https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Carolyn_Lawrence_Dill_GOMAP_Maize_MaizeGDB_B73_NAM_5.0_October_2022_v2.r1/3_final-result.
Input and output example files are in the data directory.
The script take a csv file as input and provides a txt file as output.

Usage: python zmays_edit_go_output.py  go_infile go_outiflei
'''
go_infile = sys.argv[1] 
go_outfile = sys.argv[2] # Give a output file name

with open(go_infile, 'r') as f1:
    df = pd.read_csv(f1)
    result_df = df.groupby('gene_id')['terms'].apply(lambda x: ';'.join(x)).reset_index()

with open(go_outfile, 'w') as f2:
    result_df.to_csv(f2, index=False, sep='\t')

print("Done")

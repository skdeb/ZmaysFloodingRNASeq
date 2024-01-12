import sys

'''
Editing the GO annotations format to use with GOATOOLS.
GO-term annotation data of Z. nicaraguensis, T. dactyloides, and V. cuspidata
was extracted from trinotate report file and then a trinity utility script
was used to get the GO-terms. The output file of GO-terms were than edited
for input to GOATOOLS format. The example of input and output file formats provided in the data directory.
The script take a txt file as input and provides a txt file as output.

Usage: python edit_trinotate_go_output.py go_infile go_outifle
'''
go_infile = sys.argv[1]
go_outfile = sys.argv[2] # Give a output file name

with open(go_infile, 'r') as f1:
    lines = f1.readlines()

processed_data = []

for line in lines:
    columns = line.strip().split('\t')
    if len(columns) == 2:
        protein_accession = columns[0]
        gene_ontology_terms = columns[1].replace(',', ';').replace(' ', '')
        processed_data.append(f"{protein_accession}\t{gene_ontology_terms}")

with open(go_outfile, 'w') as f2:
    for item in processed_data:
        f2.write("%s\n" % item)

print("Done")

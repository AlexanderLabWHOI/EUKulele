#!/usr/bin/python3

"""
USAGE: 

    python create-protein-table.py --infile_peptide input.fa --infile_taxonomy tax.txt --outfile_json protein-species-map.json --output output_table.txt [--delim delimiter --column column] 

Generates table useable for taxonomic placment with EUKulele. If no delimiter is provided, default is '/'. If no column header is provided, default is SOURCE_ID. 

python EUKulele/scripts/create_protein_table.py --infile_peptide EUKulele/tests/aux_data/mmetsp/reference-pep-trunc.pep.faa --infile_taxonomy  EUKulele/tests/aux_data/mmetsp/taxonomy-table.txt --outfile_json EUKulele/tests/aux_data/mmetsp/protein-map.json --output EUKulele/tests/aux_data/mmetsp/tax-table.txt --delim "/" --strain_col_id strain_name --taxonomy_col_id taxonomy --column SOURCE_ID
"""

from Bio import SeqIO
import sys
import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--infile_peptide', type = list, nargs='+', required=True) # this should be given as a list of the input files
parser.add_argument('--infile_taxonomy',  default='') # the original taxonomy file
parser.add_argument('--outfile_json', default = 'protein-species-map.json') # the protein json file to be written
parser.add_argument('--output', default = 'output_table.txt') # the output taxonomy table (formatted) to be written
parser.add_argument('--delim',  type=str, default = '/')
parser.add_argument('--strain_col_id',  type=str, default = 'strain_name') # the column which indicates the name of the strain in the taxonomy file
parser.add_argument('--taxonomy_col_id',  type=str, default = 'taxonomy') # the column which indicates the taxonomy of the strain in the taxonomy file
parser.add_argument('--column', type=str, default = 'SOURCE_ID') # can be numeric, zero-indexed, if it's a delimited part of header
# set to true if there is a column called "taxonomy" that we wish to split
parser.add_argument('--reformat_tax', dest='reformat', default=False, action='store_true') 
parser.add_argument('--euk-prot', dest='eukprot', default=False, action='store_true') # eukprot's taxonomy is too unique

args = parser.parse_args()
# if the input is a folder, you also need to add the first token in the underscore-separated name to a dictionary for that
# one - ultimately we want to know which MMETSP or whatever it came from
odict = {}
for curr_pepfile in list(args.infile_peptide):
    print("".join(curr_pepfile))
    pepfile = "".join(curr_pepfile)
    for record in SeqIO.parse(pepfile, "fasta"):
        header = record.description
        rid = record.id.replace(".","N") #record.id.split(".")[0] #record.id.replace(".","N")
        header = str(header).replace(args.delim, "hello")
        hlist = header.split("hello")
        if len(args.infile_peptide) > 1: # if there is a list of files, use the filename as the ID
            sid = pepfile.split("/")[-1].split("_")[0]
            odict[rid] = sid
        elif args.column.isdigit():
            sid = hlist[int(args.column)]  
            odict[rid] = sid
        else:
            for h in hlist:
                if args.column in h: #h.startswith(args.column):
                    sid = h.split('=')[1].strip()    
                    odict[rid] = sid
                    break
tax_file = pd.read_csv(args.infile_taxonomy, sep = "\t", encoding='latin-1')

if args.reformat:
    colnames_tax = ["Source_ID","Supergroup","Division","Class","Order","Family","Genus","Species"]
    tax_out = pd.DataFrame(columns=colnames_tax)
    for i in range(0,len(tax_file.index)):
        if not args.eukprot:
            curr_row = [tax_file[args.strain_col_id][i]] + tax_file[args.taxonomy_col_id][i].split(";")
            if len(curr_row) < (len(colnames_tax)):
                curr_row = curr_row + [""] * ((len(colnames_tax) + 1) - len(curr_row))
            elif len(curr_row) > (len(colnames_tax)):
                curr_row = curr_row[0:7] + [curr_row[8]]
            add_series = pd.Series(curr_row, index = colnames_tax)
            tax_out = tax_out.append(add_series, ignore_index=True)
        else:
            curr_row = [tax_file[args.strain_col_id][i]] + [""] * 7
            full_taxonomy = tax_file[args.taxonomy_col_id][i].split(";")
            genus_and_species = tax_file["Name_to_Use"][i].replace("_", " ")
            curr_row[1] = full_taxonomy[0] # this is generally the "supergroup" the way we have used it thus far.
            curr_row[2] = tax_file["Supergroup_UniEuk"][i] # I believe this is closest to the division
            curr_row[4] = tax_file["Taxogroup_UniEuk"][i] # I believe this is closest to the order
            curr_row[6] = tax_file["Genus_UniEuk"][i]
            curr_row[7] = genus_and_species
            print(curr_row)
            add_series = pd.Series(curr_row, index = colnames_tax)
            tax_out = tax_out.append(add_series, ignore_index=True)
            
    tax_file = tax_out
tax_file.to_csv(args.output,sep="\t")
with open(args.outfile_json, 'w') as f:
    json.dump(odict, f)

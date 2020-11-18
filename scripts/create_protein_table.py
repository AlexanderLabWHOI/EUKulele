#!/usr/bin/env python3

"""
USAGE:

    python create-protein-table.py --infile_peptide input.fa --infile_taxonomy
    tax.txt --outfile_json protein-species-map.json --output output_table.txt
    [--delim delimiter --column column]

    Generates table useable for taxonomic placment with EUKulele.
    If no delimiter is provided, default is '/'. If no column header is provided,
    default is SOURCE_ID.

    python EUKulele/scripts/create_protein_table.py --infile_peptide
    EUKulele/tests/aux_data/mmetsp/reference-pep-trunc.pep.faa --infile_taxonomy
    EUKulele/tests/aux_data/mmetsp/taxonomy-table.txt --outfile_json
    EUKulele/tests/aux_data/mmetsp/protein-map.json --output
    EUKulele/tests/aux_data/mmetsp/tax-table.txt --delim "/"
    --strain_col_id strain_name --taxonomy_col_id taxonomy --column SOURCE_ID

    create_protein_table.py --infile_peptide reference.pep.fa --infile_taxonomy
    taxonomy-table.txt --outfile_json prot-map.json --output tax-table.txt
    --delim "/" --col_source_id strain_name --taxonomy_col_id taxonomy --column 2

    python
    /vortexfs1/omics/alexander/akrinos/remodeling/EUKulele/scripts/create_protein_table.py
    --infile_peptide
    /vortexfs1/omics/alexander/akrinos/EUKulele-Reference/phylodb_db/reference.pep.fa
    --infile_taxonomy
    /vortexfs1/omics/alexander/akrinos/EUKulele-Reference/phylodb_db/tax-table.txt
    --outfile_json
    /vortexfs1/omics/alexander/akrinos/EUKulele-Reference/phylodb_db/prot-map.json
    --output
    /vortexfs1/omics/alexander/akrinos/EUKulele-Reference/phylodb_db/taxonomy_table.txt
    --delim "\t" --strain_col_id strain_name --taxonomy_col_id taxonomy --column 2

"""

import os
import argparse
import json
from Bio import SeqIO
import pandas as pd

def createProteinTable(args=None):
    '''
    Main function; intended to parse and create required files
    for EUKulele database creation.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile_peptide', type = list, nargs='+', required=True)
    # this should be given as a list of the input files
    parser.add_argument('--infile_taxonomy',  default='')
    # the original taxonomy file
    parser.add_argument('--outfile_json', default = 'prot-map.json')
    # the protein json file to be written
    parser.add_argument('--output', default = 'tax-table.txt')
    # the output taxonomy table (formatted) to be written
    parser.add_argument('--delim',  type=str, default = '/')
    parser.add_argument('--col_source_id',  type=str, default = 'Source_ID')
    # the column which indicates the name of the strain in the taxonomy file
    parser.add_argument('--taxonomy_col_id',  type=str, default = 'taxonomy')
    # the column which indicates the taxonomy of the strain in the taxonomy file
    parser.add_argument('--column', type=str, default = 'SOURCE_ID')
    # can be numeric, zero-indexed, if it's a delimited part of header
    # set to true if there is a column called "taxonomy" that we wish to split
    parser.add_argument('--reformat_tax', dest='reformat', default=False, action='store_true')
    parser.add_argument('--euk-prot', dest='eukprot',
                        default=False, action='store_true') # eukprot's taxonomy is too unique

    if args is not None:
        args = parser.parse_args(args)
    else:
        args = parser.parse_args()

    # if the input is a folder, you also need to add the first token in the
    # underscore-separated name to a dictionary for that
    # one - ultimately we want to know which MMETSP or whatever it came from

    odict = {}
    for curr_pepfile in list(args.infile_peptide):
        pepfile = "".join(curr_pepfile)
        for record in SeqIO.parse(pepfile, "fasta"):
            header = record.description
            rid = record.id.replace(".","N") #record.id.split(".")[0] #record.id.replace(".","N")
            counter = 2
            while rid in odict:
                if "_" in rid:
                    rid = "_".join(rid.split("_")[0:-1]) + "_" + str(counter)
                else:
                    rid = rid + "_" + str(counter)
                counter = counter + 1
            if 't' in args.delim: # why is this tab thing not working otherwise?? even the equality
                tester = "".join(list(str(header))).replace('\t', '    ')
                hlist = tester.split("    ")
            else:
                header = str(header).replace(args.delim, "hello")
                hlist = header.split("hello")
            if len(args.infile_peptide) > 1:
                # if there is a list of files, use the filename as the ID
                sid = pepfile.split("/")[-1].split("_")[0]
                odict[rid] = sid
            elif args.column.isdigit():
                sid = hlist[int(args.column)]
                odict[rid] = sid
            else:
                for h_curr in hlist:
                    if args.column in h_curr: #h.startswith(args.column):
                        sid = h_curr.split('=')[1].strip()
                        odict[rid] = sid
                        break

        print("Modifying...",pepfile,flush=True)
        os.system("cut -f 1 " + str(pepfile) + " > " + str(pepfile) + ".tester.pep.fa")
        os.system("perl -i -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' " + \
                  str(pepfile) + ".tester.pep.fa")
        os.system("mv " + str(pepfile) + ".tester.pep.fa " + str(pepfile))

    tax_file = pd.read_csv(args.infile_taxonomy, sep = "\t", encoding='latin-1')

    if args.reformat:
        colnames_tax = ["Source_ID","Supergroup","Division","Class",
                        "Order","Family","Genus","Species"]
        tax_out = pd.DataFrame(columns=colnames_tax)
        for i in range(0,len(tax_file.index)):
            if not args.eukprot:
                curr_row = [tax_file[args.col_source_id][i]] + \
                           tax_file[args.taxonomy_col_id][i].split(";")
                if len(curr_row) < (len(colnames_tax)):
                    curr_row = curr_row + [""] * ((len(colnames_tax) + 1) - len(curr_row))
                elif len(curr_row) > (len(colnames_tax)):
                    curr_row = curr_row[0:7] + [curr_row[8]]
                add_series = pd.Series(curr_row, index = colnames_tax)
                tax_out = tax_out.append(add_series, ignore_index=True)
            else:
                curr_row = [tax_file[args.col_source_id][i]] + [""] * 7
                full_taxonomy = tax_file[args.taxonomy_col_id][i].split(";")
                genus_and_species = tax_file["Name_to_Use"][i].replace("_", " ")
                curr_row[1] = full_taxonomy[0]
                # this is generally the "supergroup" the way we have used it thus far.
                curr_row[2] = tax_file["Supergroup_UniEuk"][i]
                # this is closest to the division
                curr_row[4] = tax_file["Taxogroup_UniEuk"][i]
                # this is closest to the order
                curr_row[6] = tax_file["Genus_UniEuk"][i]
                curr_row[7] = genus_and_species
                add_series = pd.Series(curr_row, index = colnames_tax)
                tax_out = tax_out.append(add_series, ignore_index=True)

        tax_file = tax_out
    tax_file.to_csv(args.output,sep="\t")
    with open(args.outfile_json, 'w') as file_out:
        json.dump(odict, file_out)

    return 0

if __name__ == "__main__":
    createProteinTable()

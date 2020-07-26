import os
import sys
import yaml
import subprocess
import pkgutil

import EUKulele
from EUKulele.manage_steps import createAlignmentDatabase

from scripts.create_protein_table import createProteinTable

def downloadDatabase(database_name, alignment_choice):
    """
    Automatically downloads a peptide database for use with EUKulele and stores the name of
    the resulting FASTA file and taxonomy table.
    """
    
    print("Automatically downloading database " + database_name + ". If you intended to " + \
          "use an existing database folder, be sure a reference FASTA, protein map, and taxonomy table " + \
          "are provided. Check the documentation for details.")
    create_protein_table_args = []
    if database_name == "mmetsp":
        column_id = "SOURCE_ID"
        delimiter = "/"
    elif (database_name == "eukprot"):
        column_id = 0
        delimiter = "\t"
        create_protein_table_args.append("--euk-prot")
    elif (database_name == "phylodb"):
        column_id = 0
        delimiter = "\t"
        create_protein_table_args.append("--reformat_tax")
    else:
        print("Specified reference database, " + database_name + " is not supported.")
        sys.exit(1)
        
    #config_file = pkgutil.get_data(__name__, "static/reference_url.yaml")
    config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "static/reference_url.yaml")
    with open(config_file, 'r') as configfile:
        config = yaml.safe_load(configfile)
        
    database_ref_url = config[database_name + "_ref"]
    database_tab_url = config[database_name + "_tab"]
    rc = os.system(" ".join(["source", "download_database.sh", database_name, database_ref_url, 
                          database_tab_url]))
    if rc != 0:
        print("Download of database " + database_name + " did not complete correctly.")
        sys.exit(1)
        
    fasta_name = os.path.join(database_name,"reference.pep.fa")
    orig_tax_name = os.path.join(database_name,"taxonomy-table.txt")
    
    tax_table = os.path.join(database_name,"tax-table.txt")
    protein_json = os.path.join(database_name,"tax-table.txt")
    
    create_protein_table_args.extend(["--infile_peptide",fasta_name,"--infile_taxonomy",
                                      orig_tax_name,"--output",tax_table,"--outfile_json",
                                      protein_json,"--delim",str(delimiter),"--strain_col_id",
                                      "strain_name","--taxonomy_col_id", "taxonomy","--column",
                                      str(column_id)])
    
    ## Run function to create protein table file from scripts/create_protein_table.py ##
    rc1 = createProteinTable(create_protein_table_args)
    if rc1 != 0:
        print("Taxonomy table and protein JSON file creation step did not complete successfully.")
        sys.exit(1)
        
    rc2 = createAlignmentDatabase(fasta_name, database_name)
    if rc2 != 0:
        print("Alignment database for " + alignment_choice + " did not complete successfully;" +
              "check log for details.")
        sys.exit(1)
    
    return fasta_name, tax_table, protein_json
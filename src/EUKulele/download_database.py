'''
Automates database download process.
'''

import os
import sys
import subprocess
import yaml
import datetime

import EUKulele
from EUKulele.manage_steps import createAlignmentDatabase

from scripts.create_protein_table import createProteinTable

def downloadDatabase(database_name, alignment_choice, output_dir, reference_dir = "."):
    """
    Automatically downloads a peptide database for use with EUKulele and
    stores the name of the resulting FASTA file and taxonomy table.
    """

    print("Automatically downloading database " + database_name,
          ". If you intended to " + \
          "use an existing database folder, be sure a reference FASTA,",\
          "protein map, and taxonomy table " + \
          "are provided. Check the documentation for details.")
    create_protein_table_args = []    
    
    if (database_name == "mmetsp") | (database_name == "marmmetsp"):
        column_id = "SOURCE_ID"
        delimiter = "/"
        sourceID = "Source_ID"
    elif database_name == "eukprot":
        column_id = 0
        delimiter = "\t"
        create_protein_table_args.append("--euk-prot")
        sourceID = "EukProt_ID"
    elif database_name == "phylodb":
        column_id = 2
        delimiter = "\t"
        sourceID = "strain_name"
        create_protein_table_args.append("--reformat_tax")
    elif database_name == "eukzoo":
        column_id = 0
        delimiter = "_"
        sourceID = "Source_ID"
    else:
        print("Specified reference database, " + database_name + " is not supported.")
        sys.exit(1)

    #config_file = pkgutil.get_data(__name__, "static/reference_url.yaml")
    config_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               "static/reference_url.yaml")
    with open(config_file, 'r') as configfile:
        config = yaml.safe_load(configfile)

    database_ref_url = config[database_name + "_ref"]
    database_tab_url = config[database_name + "_tab"]
    rc = os.system(" ".join(["download_database.sh", database_name,
                             database_ref_url,database_tab_url,
                             reference_dir]))
    if rc != 0:
        print("Download of database " + database_name +\
              " did not complete correctly.")
        sys.exit(1)
        
    e = datetime.datetime.now()
    f = open(os.path.join(output_dir, "README_DB.txt"), "w")
    f.write("EUKulele was run at the following time, using the following file locations:")
    f.write("Time was " + str(e) + " and database " + str(database_name) + " was downloaded.")
    f.write("The reference URL was: " + str(database_ref_url) + \
            " and the URL of the taxonomy table was: " + str(database_tab_url))
    f.close()

    fasta_name = "reference.pep.fa" #os.path.join(database_name,"reference.pep.fa")
    orig_tax_name = os.path.join(reference_dir, database_name,"taxonomy-table.txt")

    tax_table = os.path.join(reference_dir, database_name,"tax-table.txt")
    protein_json = os.path.join(reference_dir, database_name,"prot-map.json")

    create_protein_table_args.extend(["--infile_peptide",
                                      os.path.join(reference_dir,database_name,fasta_name),
                                      "--infile_taxonomy",orig_tax_name,
                                      "--output",tax_table,"--outfile_json",
                                      protein_json,"--delim",str(delimiter),
                                      "--col_source_id",sourceID,"--taxonomy_col_id",
                                      "taxonomy","--column",
                                      str(column_id)])

    ## Run function to create protein table file from scripts/create_protein_table.py ##
    createProteinTable_log = open(os.path.join(output_dir,"log",
                                               "proteintab.log"), "w+")
    createProteinTable_err = open(os.path.join(output_dir,"log",
                                               "proteintab.err"), "w+")
    sys.stdout = createProteinTable_log
    sys.stderr = createProteinTable_err
    rc1 = createProteinTable(create_protein_table_args)
    if rc1 != 0:
        print("Taxonomy table and protein JSON file creation",\
              "step did not complete successfully.")
        sys.exit(1)

    rc2 = createAlignmentDatabase(fasta_name.split("/")[-1], True, output_dir,
                                  alignment_choice,
                                  os.path.join(reference_dir,database_name))
    if rc2 != 0:
        print("Alignment database for " + alignment_choice + \
              " did not initially complete successfully; " +
              "check log (proteintab.err) for details. Trying again...")
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    return fasta_name, tax_table, protein_json

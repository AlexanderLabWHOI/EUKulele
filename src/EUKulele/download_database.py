import os
import sys
import yaml
import subprocess
import pkgutil

def createAlignmentDatabase(REF_FASTA, DATABASE_DIR, alignment_choice):
    output_log = "alignment_out.log"
    error_log = "alignment_err.log"
    if alignment_choice == "diamond":
        align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if (not os.path.isfile(align_db)) | (RERUN_RULES):
            ## DIAMOND database creation ##
            db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa'))
            rc = os.system("diamond makedb --in " + concatenated_file + " --db " + db + " 1> " + output_log + " 2> " + error_log)
        else:
            print("Diamond database file already created; will not re-create database.", flush = True)
    else:
        db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
        db_type = "prot"
        blast_version = 5
        rc = os.system("makeblastdb -in " + concatenated_file + " -parse_seqids -blastdb_version " + str(blast_version) + " -title " + database + " -dbtype " + db_type + " -out " + db)
        
    return rc

# Download the three supported eukaryote databases automatically and store the name of the resulting FASTA file
# and taxonomy table.
def downloadDatabase(database_name, alignment_choice):
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
        
    config_file = pkgutil.get_data(__name__, "src/EUKulele/static/reference_url.yaml")
    with open(config_file, 'r') as configfile:
        config = yaml.safe_load(configfile)
        
    database_ref_url = config[database_name + "_ref"]
    database_tab_url = config[database_name + "_tab"]
    p = subprocess.Popen(["source", "download_database.sh", database_name, database_ref_url, database_tab_url])
    p.wait()
    if p.returncode != 0:
        print("Download of database " + database_name + " did not complete correctly.")
        sys.exit(1)
        
    fasta_name = os.path.join(database_name,"reference.pep.fa")
    orig_tax_name = os.path.join(database_name,"taxonomy-table.txt")
    
    tax_table = os.path.join(database_name,"tax-table.txt")
    protein_json = os.path.join(database_name,"tax-table.txt")
    
    create_protein_table_args.extend(["--infile_peptide",fasta_name,"--infile_taxonomy",orig_tax_name,"--output",tax_table,"--outfile_json",protein_json,"--delim",str(delimiter),"--strain_col_id","strain_name","--taxonomy_col_id", "taxonomy","--column",str(column_id)])
    
    ## Run function to create protein table file from scripts/create_protein_table.py
    rc1 = createProteinTable(create_protein_table_args)
    if rc1 != 0:
        print("Taxonomy table and protein JSON file creation step did not complete successfully.")
        sys.exit(1)
        
    rc2 = createAlignmentDatabase(fasta_name, database_name)
    if rc2 != 0:
        print("Alignment database for " + alignment_choice + " did not complete successfully; check log for details.")
        sys.exit(1)
    
    return fasta_name, tax_table, protein_json
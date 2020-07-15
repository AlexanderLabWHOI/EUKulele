import pandas as pd
import os
import numpy as np
import math
import sys
import yaml
import argparse

def eukulele(config="", string_arguments=""):
    if (config == "") | (not os.path.isfile(config)):
        print("Running EUKulele with command line arguments, as no valid configuration file was provided.")
        os.system("python EUKulele-main.py " + str(string_arguments))
    else:
        print("Running EUKulele with entries from the provided configuration file.")
        args = parseConfig(config)
        os.system("python EUKulele-main.py " + str(args)) 

parser = argparse.ArgumentParser()
parser.add_argument("--config", dest = "configfile", default = "")
args = parser.parse_args()

if (args.configfile == "") | (not os.path.isfile(args.configfile)):
    arglist = sys.argv[1:]
    print("Running EUKulele with command line arguments, as no valid configuration file was provided.")
    os.system("python EUKulele-main.py " + str(arglist))
    sys.exit(1)

def parseConfig(configfile):
    with open(configfile, 'r') as configfile:
        config = yaml.safe_load(configfile)

    args = ""

    ## CHECK THAT ALL OPTIONS WITHOUT DEFAULTS EXIST IN CONFIGURATION FILE ##
    required_entries = ["mets_or_mags","reference","samples","cutoff"]
    for r in required_entries:
        if r not in config:
            print("You did not include required entry: " + str(r) + " in the configuration file.")
            sys.exit(1)

    ## BASIC OPTIONS ##
    args = args + " --mets_or_mags " + str(config["mets_or_mags"])
    args = args + " --reference_dir " + str(config["reference"])
    args = args + " --sample_dir " + str(config["samples"])
    if "ref_fasta" in config: # otherwise, this will default to reference.pep.fa! This works for auto-downloaded databases.
        args = args + " --ref_fasta " + str(config["ref_fasta"])
    if "output" in config:
        args = args + " --out_dir " + str(config["output"])
    if "database" in config:
        args = args + " --database " + str(config["database"])
    if "nucleotide_extension" in config:
        args = args + " --nucleotide_extension " + str(config["nucleotide_extension"])
    if "protein_extension" in config:
        args = args + " --protein_extension " + str(config["protein_extension"])
    if "download_reference" in config:
        args = args + " --download_reference"
    if "scratch" in config:
        args = args + " --scratch " + str(config["scratch"])
    if "ref_fasta_ext" in config:
        args = args + " --ref_fasta_ext " + str(config["ref_fasta_ext"])
    if "force_rerun" in config:
        if config["force_rerun"] == 1:
            args = args + " --force_rerun"

    ## SALMON OPTIONS ##
    if "use_salmon_counts" in config:
        args = args + " --use_salmon_counts " + str(config["use_salmon_counts"])
        if config["use_salmon_counts"] == 1:
            if "salmon_dir" in config:
                args = args + " --salmon_dir " + config["salmon_dir"]
            else:
                print("You need to include a salmon directory if you wish to process salmon counts.")
                sys.exit(1)
            if "names_to_reads" in config:
                args = args + " --names_to_reads " + config["names_to_reads"]

    ## TRANSDECODER AND COMPUTATIONS OPTIONS ##
    if "transdecoder_orfsize" in config:
        args = args + " --transdecoder_orfsize " + str(config["transdecoder_orfsize"])
    if "CPUs" in config:
        args = args + " --CPUs " + str(config["CPUs"])
    if "choose_parallel" in config:
        if (config["choose_parallel"] == "parallel") | (config["choose_parallel"] == 1) | (config["choose_parallel"] == "True") | (config["choose_parallel"] == True):
            args = args + " -p "

    ## ALIGNMENT AND BUSCO OPTIONS ##
    if "alignment_choice" in config: 
        args = args + " --alignment_choice " + str(config["alignment_choice"])
    if "cutoff_file" in config:    
        args = args + " --cutoff_file " + config["cutoff_file"]
    if "filter_metric" in config:
        args = args + " --filter_metric " + config["filter_metric"]
    if "consensus_cutoff" in config:
        args = args + " --consensus_cutoff " + str(config["consensus_cutoff"])
    if "busco_file" in config:
        args = args + " --busco_file " + str(config["busco_file"])
    if "individual_or_summary" in config:
        args = args + " --individual_or_summary " + str(config["individual_or_summary"])
    if "organisms" in config:
        args = args + " --organisms " + str(config["organisms"])
    if "taxonomy_organisms" in config:
        args = args + " --taxonomy_organisms " + str(config["taxonomy_organisms"])
    if "busco_threshold" in config:
        args = args + " --busco_threshold " + str(config["busco_threshold"])

    ## AUTOMATIC TAXONOMY TABLE CREATION ##
    if ("create_tax_table" in config) & (config["create_tax_table"] == 1):
        args = args + " --create_tax_table"
        if ("original_tax_table" not in config):
            print("You have specified that the program should create a formatted tax table, but did not provide an original taxonomy file.")
            sys.exit(1)
        original_tax_table = config["original_tax_table"]
        args = args + " --original_tax_table " + str(original_tax_table)

    if ("strain_col_id" in config):
        strain_col_id = config["strain_col_id"]
        args = args + " --strain_col_id " + str(strain_col_id)
    if ("taxonomy_col_id" in config):
        taxonomy_col_id = config["taxonomy_col_id"]
        args = args + " --taxonomy_col_id " + str(taxonomy_col_id)
    if ("column" in config):
        column = config["column"]
        args = args + " --column " + str(column)
    if ("reformat_tax" in config):
        column = config["reformat_tax"]
        args = args + " --reformat_tax"
    if ("delimiter" in config):
        delimiter = config["delimiter"]
        args = args + " --delimiter " + str(delimiter)
    if ("tax_table" in config): # unique, non-default name for formatted taxonomy table
        tax_table = config["tax_table"]
        args = args + " --tax_table " + str(tax_table)
    if ("protein_map" in config): # unique, non-default name for formatted protein map
        protein_map = config["protein_map"]
        args = args + " --protein_map " + str(protein_map)
        
    return args
 
euk_args = parseConfig(args.configfile)
os.system("python EUKulele-main.py " + str(euk_args))        
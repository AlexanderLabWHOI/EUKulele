import pandas as pd
import os
import numpy as np
import math
import sys
import yaml
import argparse
import pathlib
try:
    import EUKulele
except:
    pass

abs_path = os.path.abspath(os.path.dirname(__file__))

def eukulele(config="", string_arguments=""):
    split_args = string_arguments.split(" ")
    
    if len(split_args) > 2:
        if string_arguments.split(" ")[1] == "--config":
            stringargs = parseConfig(string_arguments.split(" ")[2])
            print("Running EUKulele with entries from the provided configuration file.")
            EUKulele.EUKulele_main.main(str(stringargs)) 
            return 0
        
    if (config == "") | (not os.path.isfile(config)):
        print("Running EUKulele with command line arguments, as no valid configuration file was provided.")
        
        EUKulele.EUKulele_main.main(str(string_arguments)) 
    else:
        print("Running EUKulele with entries from the provided configuration file.")
        args = parseConfig(config)
        EUKulele.EUKulele_main.main(str(args)) 
 
def eukulele_cl():
    sys.path.append(os.path.realpath('..'))
    sys.path.append(pathlib.Path(__file__).parent.absolute())
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", dest = "configfile", default = "")
    args = parser.parse_args()

    if (args.configfile == "") | (not os.path.isfile(args.configfile)):
        arglist = sys.argv[1:]
        print("Running EUKulele with command line arguments, as no valid configuration file was provided.")
        #os.system("python EUKulele_main.py " + str(arglist))
        EUKulele.EUKulele_main.main(str(arglist)) 
        sys.exit(0)
        
    euk_args = parseConfig(args.configfile)
    EUKulele.EUKulele_main.main(str(euk_args)) 
    #os.system("python EUKulele_main.py " + str(euk_args))   

def parseConfig(config_file):
    with open(config_file, 'r') as configfile:
        config = yaml.safe_load(configfile)

    args = ""

    ## CHECK THAT ALL OPTIONS WITHOUT DEFAULTS EXIST IN CONFIGURATION FILE ##
    required_entries = ["mets_or_mags","samples"]
    for r in required_entries:
        if r not in config:
            print("You did not include required entry: " + str(r) + " in the configuration file.")
            sys.exit(1)

    ## BASIC OPTIONS ##
    if "subroutine" in config:
        args = args + str(config["subroutine"]) + " "
    args = args + " --config_file " + str(config_file)
    args = args + " --mets_or_mags " + str(config["mets_or_mags"])
    
    ## If reference_dir is provided, databases are not downloaded.
    if "reference_dir" in config:
        args = args + " --reference_dir " + str(config["reference"])
    if "samples" in config:
        args = args + " --sample_dir " + str(config["samples"])
    if "ref_fasta" in config: # otherwise, this will default to reference.pep.fa! Set automatically if database is auto-downloaded ("download_reference", below)
        args = args + " --ref_fasta " + str(config["ref_fasta"])
    if "output" in config:
        args = args + " --out_dir " + str(config["output"])
    if "database" in config:
        args = args + " --database " + str(config["database"])
    if "nucleotide_extension" in config:
        args = args + " --nucleotide_extension " + str(config["nucleotide_extension"])
    if "protein_extension" in config:
        args = args + " --protein_extension " + str(config["protein_extension"])
    if "scratch" in config:
        args = args + " --scratch " + str(config["scratch"])
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

    ## ALIGNMENT AND BUSCO OPTIONS ##
    if "alignment_choice" in config: 
        args = args + " --alignment_choice " + str(config["alignment_choice"])
    if "cutoff" in config:    
        args = args + " --cutoff_file " + config["cutoff"]
    if "filter_metric" in config:
        args = args + " --filter_metric " + config["filter_metric"]
    if "consensus_cutoff" in config:
        args = args + " --consensus_cutoff " + str(config["consensus_cutoff"])
    if "busco_file" in config:
        args = args + " --busco_file " + str(config["busco_file"])
    if "individual_or_summary" in config:
        args = args + " --individual_or_summary " + str(config["individual_or_summary"])
    if "organisms" in config:
        args = args + " --organisms " + str(" ".join(config["organisms"]))
    if "taxonomy_organisms" in config:
        args = args + " --taxonomy_organisms " + str(" ".join(config["taxonomy_organisms"]))
    if "busco_threshold" in config:
        args = args + " --busco_threshold " + str(config["busco_threshold"])

    if ("strain_col_id" in config):
        strain_col_id = config["strain_col_id"]
        args = args + " --strain_col_id " + str(strain_col_id)
    if ("taxonomy_col_id" in config):
        taxonomy_col_id = config["taxonomy_col_id"]
        args = args + " --taxonomy_col_id " + str(taxonomy_col_id)
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
     
if __name__ == "__main__": 
    eukulele(string_arguments = " ".join(sys.argv))
'''
Sets up EUKulele imports for downstream processing.
'''

import os
import sys
import yaml

import EUKulele

abs_path = os.path.abspath(os.path.dirname(__file__))

def eukulele(config="", string_arguments=""):
    '''
    A function that can be called within the package or by the
    user to run EUKulele.
    '''

    split_args = string_arguments.split(" ")

    if len(split_args) > 1:
        if string_arguments.split(" ")[0] == "--config":
            stringargs = parseConfig(string_arguments.split(" ")[1])
            print("Running EUKulele with entries from the provided configuration file.")
            EUKulele.EUKulele_main.main(str(stringargs))
            return 0

    if (config == "") | (not os.path.isfile(config)):
        print("Running EUKulele with command line arguments, " + \
              "as no valid configuration file was provided.")

        EUKulele.EUKulele_main.main(str(string_arguments))
    else:
        print("Running EUKulele with entries from the provided configuration file.")
        args = parseConfig(config)
        EUKulele.EUKulele_main.main(str(args))

def parseConfig(config_file):
    '''
    Process a user-defined configuration file.
    '''

    with open(config_file, 'r') as configfile:
        config = yaml.safe_load(configfile)

    args = ""

    ## CHECK THAT ALL OPTIONS WITHOUT DEFAULTS EXIST IN CONFIGURATION FILE ##
    required_entries = ["mets_or_mags","samples"]
    for r_curr in required_entries:
        if r_curr not in config:
            print("You did not include required entry: " + \
                  str(r_curr) + " in the configuration file.")
            sys.exit(1)

    ## BASIC OPTIONS ##
    if "subroutine" in config:
        args = args + str(config["subroutine"]) + " "
    args = args + " --config_file " + str(config_file)
    args = args + " --mets_or_mags " + str(config["mets_or_mags"])

    ## If reference_dir is provided, databases are not downloaded.
    if "reference" in config:
        args = args + " --reference_dir " + str(config["reference"])
    if "samples" in config:
        args = args + " --sample_dir " + str(config["samples"])
    if "ref_fasta" in config:
        # otherwise, this will default to reference.pep.fa!
        #Set automatically if database is auto-downloaded
        #("download_reference", below)
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
        if config["use_salmon_counts"] == 1:
            args = args + " --use_salmon_counts"
            if "salmon_dir" in config:
                args = args + " --salmon_dir " + config["salmon_dir"]
            else:
                print("You need to include a salmon directory "+\
                      "if you wish to process salmon counts.")
                sys.exit(1)
            if "names_to_reads" in config:
                args = args + " --names_to_reads " + \
                       config["names_to_reads"]

    ## TRANSDECODER AND COMPUTATIONS OPTIONS ##
    if "transdecoder_orfsize" in config:
        args = args + " --transdecoder_orfsize " + \
               str(config["transdecoder_orfsize"])
    if "CPUs" in config:
        args = args + " --CPUs " + str(config["CPUs"])
    if "run_transdecoder" in config:
        args = args + " --run_transdecoder"

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
        args = args + " --individual_or_summary " + \
               str(config["individual_or_summary"])
    if "organisms" in config:
        args = args + " --organisms " + \
               str(" ".join(config["organisms"]))
    if "taxonomy_organisms" in config:
        args = args + " --taxonomy_organisms " + \
               str(" ".join(config["taxonomy_organisms"]))
    if "busco_threshold" in config:
        args = args + " --busco_threshold " + str(config["busco_threshold"])
    if "test" in config:
        if config["test"] == 1:
            args = args + " logs/cdhit/mega_merge_err.log"

    if "tax_table" in config: # unique, non-default name for formatted taxonomy table
        tax_table = config["tax_table"]
        args = args + " --tax_table " + str(tax_table)
    if "protein_map" in config: # unique, non-default name for formatted protein map
        protein_map = config["protein_map"]
        args = args + " --protein_map " + str(protein_map)

    return args

if __name__ == "__main__":
    eukulele(string_arguments = " ".join(sys.argv))

#!/usr/bin/env python

# Software for taxonomic identification of eukaryotes.

import pandas as pd
import numpy as np
import os
import sys
import yaml
import chardet
import argparse
import multiprocessing
import subprocess
import shutil
import glob
from joblib import Parallel, delayed
import pkgutil

import tax_placement
from tax_placement import *

import visualize_results
from visualize_results import *

import EUKulele
from EUKulele.download_database import downloadDatabase
from EUKulele.manage_steps import manageEukulele
from EUKulele.busco_runner import readBuscoFile

import scripts as HelperScripts
from scripts.mag_stats import magStats
from scripts.names_to_reads import namesToReads
from scripts.query_busco import queryBusco

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"


def assign_taxonomy(sample_name, OUTPUTDIR, mets_or_mags):
    taxfile = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "-estimated-taxonomy.out")
    levels_directory = os.path.join(OUTPUTDIR, mets_or_mags, "levels")
    max_dir = os.path.join(OUTPUTDIR, mets_or_mags)
    error_log = os.path.join("log", "tax_assign_" + sample_name + ".err")
    out_log = os.path.join("log", "tax_assign_" + sample_name + ".out")
    
    sys.stdout = open(out_log, "w")
    sys.stderr = open(error_log, "w")
    rc = magStats(["--estimated-taxonomy-file",taxfile,"--out-prefix",sample_name,"--outdir",levels_directory,"--max-out-dir",max_dir])
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return rc
    
def run_busco(sample_name, outputdir, busco_db, mets_or_mags, PEP_EXT, NT_EXT, CPUS, OUTPUTDIR, SAMPLE_DIR):
    if mets_or_mags == "mets":
        fastaname = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
    else:
        fastaname = os.path.join(SAMPLE_DIR, sample_name + "." + PEP_EXT)
        
    rc2 = 0
    busco_run_log = os.path.join("log","busco_run.out")
    busco_run_err = os.path.join("log","busco_run.err")
    busco_config_log = os.path.join("log","busco_config.out")
    busco_config_err = os.path.join("log","busco_config.err")
    
    if not os.path.isdir(os.path.join("busco_downloads","lineages","eukaryota_odb10")):
        rc2 = os.system(" ".join(["configure_busco.sh", str(sample_name), str(outputdir), outputdir + "/config_" + sample_name + ".ini", fastaname, str(CPUS), busco_db]) + " 1> " + busco_config_log + " 2> " + busco_config_err)
    rc1 = os.system(" ".join(["run_busco.sh", str(sample_name), str(outputdir), outputdir + "/config_" + sample_name + ".ini", fastaname, str(CPUS), busco_db]) + " 1> " + busco_run_log + " 2> " + busco_run_err)
    return rc1 + rc2

def main(args_in):
    parser = argparse.ArgumentParser(
        description='Thanks for using EUKulele! EUKulele is a standalone taxonomic annotation software.\n'
                    'EUKulele is designed primarily for marine microbial eukaryotes. Check the README '
                    'for further information.',
        usage='eukulele [subroutine] --mets_or_mags [dataset_type] --reference_dir [reference_database_location] ' + 
              '--sample_dir [sample_directory] [all other options]')
    
    parser.add_argument('subroutine', metavar="subroutine", nargs='?', type=str, default="all", 
                        choices = ["","all","setup","alignment","busco"], help='Choice of subroutine to run.')
    parser.add_argument('--mets_or_mags', required = True) 
    parser.add_argument('--n_ext', '--nucleotide_extension', dest = "nucleotide_extension", default = ".fasta") 
    parser.add_argument('--p_ext', '--protein_extension', dest = "protein_extension", default = ".faa") 
    parser.add_argument('f', '--force_rerun', action='store_true', default=False)
    parser.add_argument('--scratch', default = '../scratch', 
                        help = "The scratch location to store intermediate files.")
    parser.add_argument('--config_file', default = '')

    ## SALMON OPTIONS ##
    parser.add_argument('--use_salmon_counts', type = int, default = 0)
    parser.add_argument('--salmon_dir', 
                        help = "Salmon directory is required if use_salmon_counts is true.")
    parser.add_argument('--names_to_reads',default=0, help = "A file to be created or used if it exists " +
                        "that relates transcript names to salmon counts from the salmon directory.")

    ## WHERE FILES ARE LOCATED ##
    parser.add_argument('--database', default="mmetsp", 
                        help = "The name of the database to be used to assess the reads.")
    parser.add_argument('-o','--out_dir', dest = "out_dir", default = "output", 
                        help = "Folder where the output will be written.")
    parser.add_argument('--sample_dir', required = True, 
                        help = "Folder where the input data is located (the protein or peptide files to be assessed).")
    
    ## ONLY SPECIFY THESE ARGUMENTS IF YOU HAVE ALREADY PROVIDED AND FORMATTED YOUR OWN DATABASE ##
    parser.add_argument('--reference_dir', default="", 
                        help = "Folder containing the reference files for the chosen database.")
    parser.add_argument('--ref_fasta', default = "reference.pep.fa", 
                        help = "Either a file in the reference directory where the fasta file for the database " + 
                               "is located, or a directory containing multiple fasta files that " + 
                               "constitute the database.")
  
    parser.add_argument('--tax_table', default = "tax-table.txt")
    parser.add_argument('--protein_map', default = "protein-map.json")
    
    ## ALIGNMENT OPTIONS ##
    parser.add_argument('--alignment_choice', default = "diamond", choices = ["diamond", "blast"])

    ## OPTIONS FOR CHECKING BUSCO COMPLETENESS FOR TAXONOMY ##
    parser.add_argument('--busco_file', default = "", type = str, 
                        help = "If specified, the following two arguments ('--organisms' and " + 
                        "'--taxonomy_organisms' are overwritten by the two columns of this tab-separated file.")
    parser.add_argument('--individual_or_summary', '-i', default="summary", choices=["summary","individual"], 
                        help = "These arguments are used if 'individual' is specified.") 
    parser.add_argument('--organisms', default = "", nargs = "+", 
                        help = "List of organisms to check BUSCO completeness on.")
    parser.add_argument('--taxonomy_organisms', default = "", nargs = "+", 
                        help = "Taxonomic level of organisms specified in organisms tag.")

    ## OTHER USER CHOICES ## 
    cutoff_file = pkgutil.get_data(__name__, "src/EUKulele/static/tax-cutoffs.yaml")
    parser.add_argument('--cutoff_file', default = cutoff_file)
    parser.add_argument('--filter_metric', default = "evalue", choices = ["pid", "evalue", "bitscore"])
    parser.add_argument('--consensus_cutoff', default = 0.75, type = float)
    parser.add_argument('--transdecoder_orfsize', default = 100, type = int)

    parser.add_argument('--CPUs', default=multiprocessing.cpu_count())
    parser.add_argument('--busco_threshold', default=50)
               
    args = parser.parse_args(list(filter(None, args_in.split(" "))))
    if args.subroutine == "":
        args.subroutine = "all"
    
    ## VARIABLES ##
    CONSENSUS_CUTOFF = args.consensus_cutoff
    REFERENCE_DIR = args.reference_dir
    OUTPUTDIR = args.out_dir
    SAMPLE_DIR = args.sample_dir
    REF_FASTA = args.ref_fasta

    TAX_TAB = os.path.join(REFERENCE_DIR, args.tax_table)
    PROT_TAB = os.path.join(REFERENCE_DIR, args.protein_map)
    ALIGNMENT_CHOICE = args.alignment_choice
    OUTPUT_EXTENSION = "txt"
    DBEXTENSION = ""
    TRANSDECODERORFSIZE=args.transdecoder_orfsize
    if ALIGNMENT_CHOICE == "diamond":
        OUTPUT_EXTENSION = "out"
        DBEXTENSION = ".dmnd"
    NT_EXT = args.nucleotide_extension.strip('.')
    PEP_EXT = args.protein_extension.strip('.')
    mets_or_mags = args.mets_or_mags.lower()
    
    if (mets_or_mags != "mets") & (mets_or_mags != "mags"):
        print("Only METs or MAGs are supported as input data types. Please update the 'mets_or_mags' flag accordingly.")
        sys.exit(1)

    USE_SALMON_COUNTS = args.use_salmon_counts
    SALMON_DIR = args.salmon_dir
    NAMES_TO_READS = os.path.join(REFERENCE_DIR, str(args.names_to_reads))
    CPUS = args.CPUs                         

    ORGANISMS = args.organisms 
    ORGANISMS_TAXONOMY = args.taxonomy_organisms
    BUSCO_FILE = args.busco_file
    RERUN_RULES = args.force_rerun
    
    ORGANISMS, ORGANISMS_TAXONOMY = readBuscoFile(individual_or_summary, busco_file, 
                                                  organisms, organisms_taxonomy)

    SETUP = False
    ALIGNMENT = False
    BUSCO = False
    if (args.subroutine == "all") | (args.subroutine == "setup"):
        SETUP = True
    if (args.subroutine == "all") | (args.subroutine == "alignment"):
        ALIGNMENT = True
    if (args.subroutine == "all") | (args.subroutine == "busco"):
        BUSCO = True

    ## SETUP STEPS / DOWNLOAD DEPENDENCIES ##
    manageEukulele(piece = "setup_eukulele)
    samples = manageEukulele(piece = "get_samples", mets_or_mags = mets_or_mags)

    ## Download the reference database if specified.
    if (args.reference_dir == "") | (not os.path.isdir(args.referencedir)) | \
        (not os.path.isdir(os.path.join(args.referencedir), REF_FASTA)):
        REF_FASTA, TAX_TAB, PROT_TAB = downloadDatabase(args.database.lower(), args.alignment_choice)
        REFERENCE_DIR = args.database.lower()

    if SETUP:
        manageEukulele(piece = "setup_databases", ref_fasta = REF_FASTA, rerun_rules = RERUN_RULES, 
                       alignment_choice = ALIGNMENT_CHOICE, database_dir = DATABASE_DIR)

    if ALIGNMENT:
        ## First, we need to perform TransDecoder if needed
        manageEukulele(piece = "transdecode", mets_or_mags = mets_or_mags)
        
        ## Next to align against our database of choice ##
        alignment_res = manageAlignment(alignment_choice = ALIGNMENT_CHOICE, sample_names = samples, 
                                        filter_metric = args.filter_metric, output_dir = OUTPUTDIR, 
                                        ref_fasta = REF_FASTA, mets_or_mags = mets_or_mags, database_dir = DATABASE_DIR,
                                        sample_dir = SAMPLE_DIR, rerun_rules = RERUN_RULES, 
                                        nt_ext = NT_EXT, pep_ext = PEP_EXT)

        ## Next to do salmon counts estimation (NOTE: currently only supported with a config file). ##
        if (USE_SALMON_COUNTS == 1):
            if args.config_file != "":
                rc1 = namesToReads(args.config_file)

        manageEukulele(piece = "estimate_taxonomy", output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, 
                       tax_tab = TAX_TAB, args.cutoff_file, consensus_cutoff = CONSENSUS_CUTOFF, 
                       prot_tab = PROT_TAB, use_salmon_counts = USE_SALMON_COUNTS, 
                       names_to_reads = NAMES_TO_READS, alignment_res = alignment_res, 
                       rerun_rules = RERUN_RULES)

        ## Now to visualize the taxonomy ##
        print("Performing taxonomic visualization steps...", flush=True)
        out_prefix = OUTPUTDIR.split("/")[-1]
        visualize_all_results(out_prefix, OUTPUTDIR, os.path.join(OUTPUTDIR, mets_or_mags), SAMPLE_DIR, PEP_EXT, NT_EXT, USE_SALMON_COUNTS, RERUN_RULES)

        if mets_or_mags == "mags":
            print("Performing taxonomic assignment steps...", flush=True)
            n_jobs_viz = min(multiprocessing.cpu_count(), len(samples))
            assign_res = Parallel(n_jobs=n_jobs_viz, prefer="threads")(delayed(assign_taxonomy)(samp, OUTPUTDIR, mets_or_mags) for samp in samples)
            if sum(assign_res) != 0:
                print("Taxonomic assignment of MAGs did not complete successfully. Check log files for details.")
                sys.exit(1)

    if BUSCO:
        print("Performing BUSCO steps...", flush=True)
        print("Running busco...")
        ## Run BUSCO on the full dataset ##
        busco_db = "eukaryota_odb10"
        busco_res = Parallel(n_jobs=multiprocessing.cpu_count(), prefer="threads")(delayed(run_busco)(sample_name, os.path.join(OUTPUTDIR, "busco"), busco_db, mets_or_mags, PEP_EXT, NT_EXT, CPUS, OUTPUTDIR, SAMPLE_DIR) for sample_name in samples)
        all_codes = sum(busco_res)
        if all_codes > 0:
            print("BUSCO initial run or configuration did not complete successfully. Please check the BUSCO run and configuration log files in the log/ folder.")
            sys.exit(1)

        ## Assess BUSCO completeness on the most prevalent members of the metatranscriptome at each taxonomic level ##
        if args.individual_or_summary == "individual":
            for sample_name in samples:
                busco_table = os.path.join(OUTPUTDIR, "busco", sample_name, "full_table.tsv") # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
                taxfile_stub = os.path.join(OUTPUTDIR,OUTPUTDIR.split("/")[-1]) # the prefix to specify where the taxonomy estimation output files are located

                if mets_or_mags == "mets":
                    fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
                else:
                    fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)

                query_busco_log = open(os.path.join("log","busco_query_" + sample_name + ".log"), "w+")
                query_busco_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
                sys.stdout = query_busco_log
                sys.stderr = query_busco_err
                query_args = ["--organism_group",str(" ".join(ORGANISMS)),"--taxonomic_level",str(" ".join(ORGANISMS_TAXONOMY)),"--output_dir",OUTPUTDIR,"--fasta_file",fasta,"--sample_name",sample_name,"--taxonomy_file_prefix",taxfile_stub,"--tax_table",TAX_TAB,"--busco_out",busco_table,"-i","individual"]
                rc = queryBusco(query_args)
                
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__ 
                if rc != 0:
                    print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")
        else:
            for sample_name in samples:
                busco_table = os.path.join(OUTPUTDIR, "busco", sample_name, "full_table.tsv") # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
                taxfile_stub = os.path.join(OUTPUTDIR,OUTPUTDIR.split("/")[-1]) # the prefix to specify where the taxonomy estimation output files are located

                if mets_or_mags == "mets":
                    fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
                else:
                    fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)

                query_busco_log = open(os.path.join("log","busco_query_" + sample_name + ".log"), "w+")
                query_busco_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
                sys.stdout = query_busco_log
                sys.stderr = query_busco_err
                query_args = ["--output_dir",OUTPUTDIR,"--fasta_file",fasta,"--sample_name",sample_name,"--taxonomy_file_prefix",taxfile_stub,"--tax_table",TAX_TAB,"--busco_out",busco_table,"-i","summary"]
                
                rc = queryBusco(query_args)
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__ 
                
                if rc != 0:
                    print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")

if __name__ == "__main__": 
    main(args_in = " ".join(sys.argv[1:]))
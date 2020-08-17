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

import EUKulele

from EUKulele.download_database import downloadDatabase
from EUKulele.manage_steps import manageEukulele
from EUKulele.busco_runner import readBuscoFile
from EUKulele.busco_runner import configRunBusco
from EUKulele.busco_runner import manageBuscoQuery

import scripts as HelperScripts
from scripts.names_to_reads import namesToReads

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

def main(args_in):
    parser = argparse.ArgumentParser(
        description='Thanks for using EUKulele! EUKulele is a standalone taxonomic annotation software.\n'
                    'EUKulele is designed primarily for marine microbial eukaryotes. Check the README '
                    'for further information.',
        usage='eukulele [subroutine] --mets_or_mags [dataset_type] --sample_dir [sample_directory] ' + 
              '--reference_dir [reference_database_location] [all other options]')
    
    parser.add_argument('subroutine', metavar="subroutine", nargs='?', type=str, default="all", 
                        choices = ["","all","setup","alignment","busco","coregenes"], 
                        help='Choice of subroutine to run.')
    parser.add_argument('-m', '--mets_or_mags', dest = "mets_or_mags", required = True) 
    parser.add_argument('--n_ext', '--nucleotide_extension', dest = "nucleotide_extension", default = ".fasta") 
    parser.add_argument('--p_ext', '--protein_extension', dest = "protein_extension", default = ".faa") 
    parser.add_argument('-f', '--force_rerun', action='store_true', default=False)
    parser.add_argument('--scratch', default = '../scratch', 
                        help = "The scratch location to store intermediate files.")
    parser.add_argument('--config_file', default = '')

    ## SALMON OPTIONS ##
    parser.add_argument('--use_salmon_counts', action='store_true', default=False)
    parser.add_argument('--salmon_dir', 
                        help = "Salmon directory is required if use_salmon_counts is true.")
    parser.add_argument('--names_to_reads',default=0, help = "A file to be created or used if it exists " +
                        "that relates transcript names to salmon counts from the salmon directory.")

    ## WHERE FILES ARE LOCATED ##
    parser.add_argument('-d', '--database', default="mmetsp", 
                        help = "The name of the database to be used to assess the reads.")
    parser.add_argument('-o', '--out_dir', dest = "out_dir", default = "output", 
                        help = "Folder where the output will be written.")
    parser.add_argument('-s', '--sample_dir', required = True, dest = "sample_dir",
                        help = "Folder where the input data is located (the protein or peptide files to be assessed).")
    
    ## ONLY SPECIFY THESE ARGUMENTS IF YOU HAVE ALREADY PROVIDED AND FORMATTED YOUR OWN DATABASE ##
    parser.add_argument('--reference_dir', default="", 
                        help = "Folder containing the reference files for the chosen database.")
    parser.add_argument('--ref_fasta', default = "reference.pep.fa", 
                        help = "Either a file in the reference directory where the fasta file for the database " + 
                               "is located, or a directory containing multiple fasta files that " + 
                               "constitute the database.")
  
    parser.add_argument('--tax_table', default = "tax-table.txt")
    parser.add_argument('--protein_map', default = "prot-map.json")
    
    ## ALIGNMENT OPTIONS ##
    parser.add_argument('--alignment_choice', default = "diamond", choices = ["diamond", "blast"])

    ## OPTIONS FOR CHECKING BUSCO COMPLETENESS FOR TAXONOMY ##
    parser.add_argument('--busco_file', default = "", type = str, 
                        help = "If specified, the following two arguments ('--organisms' and " + 
                        "'--taxonomy_organisms' are overwritten by the two columns of this tab-separated file.")
    parser.add_argument('--individual_or_summary', default="summary", choices=["summary","individual"], 
                        help = "These arguments are used if 'individual' is specified.") 
    parser.add_argument('-i', '--individual', dest = "individual_tag", action='store_true', default=False)
    parser.add_argument('--organisms', default = "", nargs = "+", 
                        help = "List of organisms to check BUSCO completeness on.")
    parser.add_argument('--taxonomy_organisms', default = "", nargs = "+", 
                        help = "Taxonomic level of organisms specified in organisms tag.")
    parser.add_argument('--test', action='store_true', default=False)

    ## OTHER USER CHOICES ## 
    #cutoff_file = pkgutil.get_data(__name__, "tax-cutoffs.yaml")
    cutoff_file = "tax-cutoffs.yaml"
    parser.add_argument('--cutoff_file', default = cutoff_file)
    parser.add_argument('--filter_metric', default = "evalue", choices = ["pid", "evalue", "bitscore"])
    parser.add_argument('--consensus_cutoff', default = 0.75, type = float)
    parser.add_argument('--transdecoder_orfsize', default = 100, type = int)

    parser.add_argument('--CPUs', default=multiprocessing.cpu_count())
    parser.add_argument('--busco_threshold', default=50)
    parser.add_argument('--create_fasta', action='store_true', default=False, 
                       help = "Whether to create FASTA files containing ID'd transcripts during BUSCO analysis.")
    parser.add_argument('--run_transdecoder', action='store_true', default=False,
                       help = "Whether TransDecoder should be run on metatranscriptomic samples. Otherwise, " +
                       "BLASTp is run if protein translated samples are provided, otherwise BLASTx is run " + 
                       "on nucleotide samples.")
               
    args = parser.parse_args(list(filter(None, args_in.split(" "))))
    
    ## VARIABLES ##
    TEST = args.test
    CONSENSUS_CUTOFF = args.consensus_cutoff
    REFERENCE_DIR = args.reference_dir
    OUTPUTDIR = args.out_dir
    SAMPLE_DIR = args.sample_dir
    REF_FASTA = args.ref_fasta

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
    individual_or_summary = args.individual_or_summary
    if args.individual_tag:
        individual_or_summary = "individual"
    
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
    RUN_TRANSDECODER = args.run_transdecoder
    
    ORGANISMS, ORGANISMS_TAXONOMY = readBuscoFile(individual_or_summary, BUSCO_FILE, 
                                                  ORGANISMS, ORGANISMS_TAXONOMY)

    SETUP = False
    ALIGNMENT = False
    BUSCO = False
    COREGENES = False
    if (args.subroutine == "all") | (args.subroutine == "setup"):
        SETUP = True
    if (args.subroutine == "all") | (args.subroutine == "alignment"):
        ALIGNMENT = True
    if (args.subroutine == "all") | (args.subroutine == "busco"):
        BUSCO = True
    if (args.subroutine == "all") | (args.subroutine == "coregenes"):
        COREGENES = True
        
    if TEST:
        SETUP = False
        ALIGNMENT = False
        BUSCO = False
        COREGENES = False
    else: 
        ## SETUP STEPS / DOWNLOAD DEPENDENCIES ##
        manageEukulele(piece = "setup_eukulele", output_dir = OUTPUTDIR)
        samples = manageEukulele(piece = "get_samples", mets_or_mags = mets_or_mags,
                                 sample_dir = SAMPLE_DIR, nt_ext = NT_EXT, pep_ext = PEP_EXT)

        TAX_TAB = os.path.join(REFERENCE_DIR, args.tax_table)
        PROT_TAB = os.path.join(REFERENCE_DIR, args.protein_map)

        ## Download the reference database if specified.
        if (not os.path.isfile(os.path.join(REFERENCE_DIR, REF_FASTA))) | \
           (not os.path.isfile(TAX_TAB)) | \
           (not os.path.isfile(PROT_TAB)):
            REFERENCE_DIR = args.database.lower()
            print("Specified reference directory, reference FASTA, and protein map/taxonomy table not found. " +
                  "Using database: " + REFERENCE_DIR + ".")
            os.system("mkdir -p " + REFERENCE_DIR)
            TAX_TAB = os.path.join(REFERENCE_DIR, "tax-table.txt")
            PROT_TAB = os.path.join(REFERENCE_DIR, "prot-map.json")

        if (not os.path.isfile(os.path.join(REFERENCE_DIR, REF_FASTA))) | \
           (not os.path.isfile(TAX_TAB)) | \
           (not os.path.isfile(PROT_TAB)):
            REF_FASTA, TAX_TAB, PROT_TAB = downloadDatabase(args.database.lower(), ALIGNMENT_CHOICE)
        else:
            print("Found database folder for " + REFERENCE_DIR + " in current directory; will not re-download.")

    if SETUP:
        manageEukulele(piece = "setup_databases", ref_fasta = REF_FASTA, rerun_rules = RERUN_RULES, 
                       alignment_choice = ALIGNMENT_CHOICE, database_dir = REFERENCE_DIR)

    if ALIGNMENT:
        ## First, we need to perform TransDecoder if needed
        manageEukulele(piece = "transdecode", mets_or_mags = mets_or_mags, samples = samples, output_dir = OUTPUTDIR, 
                       rerun_rules = RERUN_RULES, sample_dir = SAMPLE_DIR, transdecoder_orf_size = TRANSDECODERORFSIZE, 
                       nt_ext = NT_EXT, pep_ext = PEP_EXT, run_transdecoder = RUN_TRANSDECODER)
        
        ## Next to align against our database of choice ##
        alignment_res = manageEukulele(piece = "align_to_db", alignment_choice = ALIGNMENT_CHOICE, samples = samples, 
                                        filter_metric = args.filter_metric, output_dir = OUTPUTDIR, 
                                        ref_fasta = REF_FASTA, mets_or_mags = mets_or_mags, database_dir = REFERENCE_DIR,
                                        sample_dir = SAMPLE_DIR, rerun_rules = RERUN_RULES, 
                                        nt_ext = NT_EXT, pep_ext = PEP_EXT)

        ## Next to do salmon counts estimation. ##
        if (USE_SALMON_COUNTS == 1):
            try:
                NAMES_TO_READS = namesToReads(args.config_file)
            except:
                print("The salmon directory provided could not be converted to a salmon file. Check " +
                      "above error messages. EUKulele will continue running without counts.")
                USE_SALMON_COUNTS = 0

        manageEukulele(piece = "estimate_taxonomy", output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, 
                       tax_tab = TAX_TAB, cutoff_file = args.cutoff_file, 
                       consensus_cutoff = CONSENSUS_CUTOFF, prot_tab = PROT_TAB, use_salmon_counts = USE_SALMON_COUNTS, 
                       names_to_reads = NAMES_TO_READS, alignment_res = alignment_res, 
                       rerun_rules = RERUN_RULES, samples = samples, sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT,
                       nt_ext = NT_EXT)

        ## Now to visualize the taxonomy ##
        manageEukulele(piece = "visualize_taxonomy", output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, 
                       sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT, nt_ext = NT_EXT, 
                       use_salmon_counts = USE_SALMON_COUNTS, rerun_rules = RERUN_RULES)

        ## Next to assign taxonomy ##
        manageEukulele(piece = "assign_taxonomy", samples = samples, mets_or_mags = mets_or_mags, 
                       sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT,
                       output_dir = OUTPUTDIR)

    if BUSCO:
        configRunBusco(output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, pep_ext = PEP_EXT, 
                       nt_ext = NT_EXT, sample_dir = SAMPLE_DIR, samples = samples)

        manageBuscoQuery(output_dir = OUTPUTDIR, individual_or_summary = individual_or_summary, 
                         samples = samples, mets_or_mags = mets_or_mags, pep_ext = PEP_EXT, 
                         nt_ext = NT_EXT, sample_dir = SAMPLE_DIR, organisms = ORGANISMS, 
                         organisms_taxonomy = ORGANISMS_TAXONOMY, tax_tab = TAX_TAB, 
                         busco_threshold = args.busco_threshold)
        
    if COREGENES:
        print("Investigating core genes...")
        ## Next to align against our database of choice ##
        alignment_res = manageEukulele(piece = "core_align_to_db", alignment_choice = ALIGNMENT_CHOICE, samples = samples, 
                                        filter_metric = args.filter_metric, output_dir = OUTPUTDIR, 
                                        ref_fasta = REF_FASTA, mets_or_mags = mets_or_mags, database_dir = REFERENCE_DIR,
                                        sample_dir = SAMPLE_DIR, rerun_rules = RERUN_RULES, 
                                        nt_ext = NT_EXT, pep_ext = PEP_EXT)
        
        manageEukulele(piece = "core_estimate_taxonomy", output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, 
                       tax_tab = TAX_TAB, cutoff_file = args.cutoff_file, 
                       consensus_cutoff = CONSENSUS_CUTOFF, prot_tab = PROT_TAB, use_salmon_counts = USE_SALMON_COUNTS, 
                       names_to_reads = NAMES_TO_READS, alignment_res = alignment_res, 
                       rerun_rules = RERUN_RULES, samples = samples, sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT,
                       nt_ext = NT_EXT)

        ## Now to visualize the taxonomy ##
        manageEukulele(piece = "core_visualize_taxonomy", output_dir = OUTPUTDIR, mets_or_mags = mets_or_mags, 
                       sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT, nt_ext = NT_EXT, 
                       use_salmon_counts = USE_SALMON_COUNTS, rerun_rules = RERUN_RULES)

        ## Next to assign taxonomy ##
        manageEukulele(piece = "core_assign_taxonomy", samples = samples, mets_or_mags = mets_or_mags, 
                       sample_dir = SAMPLE_DIR, pep_ext = PEP_EXT,
                       output_dir = OUTPUTDIR)
        
        
    print("EUKulele run complete!", flush = True)
                   
if __name__ == "__main__": 
    main(args_in = " ".join(sys.argv[1:]))
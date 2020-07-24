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

abs_path = os.path.abspath(os.path.dirname(__file__))

sys.path.insert(1, os.path.join(abs_path))
sys.path.insert(1, os.path.join(abs_path, "..", "..", "scripts"))

import tax_placement
from tax_placement import *

import visualize_results
from visualize_results import *

import EUKulele
from EUKulele.setup_EUKulele import setupEukulele
from EUKulele.setup_EUKulele import setupDatabases
from EUKulele.transdecode_to_peptide import transdecodeToPeptide
from EUKulele.download_database import downloadDatabase

import scripts as HelperScripts
from scripts.create_protein_table import createProteinTable
from scripts.mag_stats import magStats
from scripts.names_to_reads import namesToReads
from scripts.query_busco import queryBusco

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

def align_to_database(alignment_choice, sample_name, filter_metric, OUTPUTDIR, REF_FASTA, mets_or_mags, DATABASE_DIR, SAMPLE_DIR, RERUN_RULES, NT_EXT, PEP_EXT):
    if alignment_choice == "diamond":
        os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "diamond"))
        diamond_out = os.path.join(OUTPUTDIR, mets_or_mags, "diamond", sample_name + ".diamond.out")
        if (os.path.isfile(diamond_out)) & (not RERUN_RULES):
            print("Diamond alignment file already detected; will not re-run step.")
            return diamond_out
        
        align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if mets_or_mags == "mets":
            fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
        else:
            fasta = os.path.join(SAMPLE_DIR, sample_name + "." + PEP_EXT)
        other = "--outfmt 6 -k 100 -e 1e-5"
        outfmt = 6
        k = 100
        e = 1e-5
        bitscore = 50
        diamond_log = os.path.join("log","diamond_align_" + sample_name + ".log")
        diamond_err = os.path.join("log","diamond_align_" + sample_name + ".err")
        if filter_metric == "bitscore":
            rc1 = os.system("diamond blastp --db " + align_db + " -q " + fasta + " -o " + diamond_out + " --outfmt " + str(outfmt) + " -k " + str(k) + " --min-score " + str(bitscore) + " 1> " + diamond_log + " 2> " + diamond_err)
        else:
            rc1 = os.system("diamond blastp --db " + align_db + " -q " + fasta + " -o " + diamond_out + " --outfmt " + str(outfmt) + " -k " + str(k) + " -e " + str(e) + " 1> " + diamond_log + " 2> " + diamond_err)
        if rc1 != 0:
            print("Diamond did not complete successfully.")
            os.system("rm -f " + diamond_out)
            return 1
        return diamond_out
    else:
        blast_out = os.path.join(OUTPUTDIR, mets_or_mags, "blast", sample_name + ".blast.txt")
        if (os.path.isfile(blast_out)) & (not RERUN_RULES):
            print("BLAST alignment file already detected; will not re-run step.")
            return blast_out
        align_db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
        if mets_or_mags == "mets":
            fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
        else:
            fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)
        outfmt = 6 # tabular output format
        e = 1e-5
        os.system("export BLASTDB=" + align_db)
        rc1 = os.system("blastp -query " + align_db + " -db " + align_db + " -out " + blast_out + " -outfmt " + str(outfmt) + " -evalue " + str(e))
        if rc1 != 0:
            print("BLAST did not complete successfully.")
            return 1
        return blast_out

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
        usage='eukulele [subroutine] --mets_or_mags [dataset_type] --reference_dir [reference_database_location] --sample_dir [sample_directory] [all other options]')
    
    parser.add_argument('subroutine', metavar="subroutine", nargs='?', type=str, default="all", choices = ["","all","setup","alignment","busco"], help='Choice of subroutine to run.')
    parser.add_argument('--mets_or_mags', required = True) 
    parser.add_argument('--n_ext', '--nucleotide_extension', dest = "nucleotide_extension", default = ".fasta") 
    parser.add_argument('--p_ext', '--protein_extension', dest = "protein_extension", default = ".faa") 
    parser.add_argument('f', '--force_rerun', action='store_true', default=False)
    parser.add_argument('--scratch', default = '../scratch') # the scratch location to store intermediate files
    parser.add_argument('--config_file', default = '')

    ## SALMON OPTIONS ##
    parser.add_argument('--use_salmon_counts', type = int, default = 0)
    parser.add_argument('--salmon_dir') # salmon directory is required if use_salmon_counts is true.
    parser.add_argument('--names_to_reads',default=0) # a file to be created or used if it exists that relates transcript names to salmon counts from the salmon directory 

    ## WHERE FILES ARE LOCATED ##
    parser.add_argument('--database', default="mmetsp") # the name of the database to be used to assess the reads
    parser.add_argument('-o','--out_dir', dest = "out_dir", default = "output") # folder where the output will be written
    parser.add_argument('--sample_dir', required = True) # folder where the input data is located (the protein or peptide files to be assessed)
    
    ## ONLY SPECIFY THESE ARGUMENTS IF YOU HAVE ALREADY PROVIDED AND FORMATTED YOUR OWN DATABASE ##
    parser.add_argument('--reference_dir', default="") # folder containing the reference files for the chosen database
    parser.add_argument('--ref_fasta', default = "reference.pep.fa") # either a file in the reference directory where the fasta file for the database is located, or a directory containing multiple fasta files that constitute the database.
    parser.add_argument('--tax_table', default = "tax-table.txt")
    parser.add_argument('--protein_map', default = "protein-map.json")
    
    ## ALIGNMENT OPTIONS ##
    parser.add_argument('--alignment_choice', default = "diamond", choices = ["diamond", "blast"])

    ## OPTIONS FOR CHECKING BUSCO COMPLETENESS FOR TAXONOMY ##
    parser.add_argument('--busco_file', default = "", type = str) # if specified, the following two arguments ("--organisms" and "--taxonomy_organisms" are overwritten by the two columns of this tab-separated file
    parser.add_argument('--individual_or_summary','-i',default="summary",choices=["summary","individual"])
    # These arguments are used if "individual" is specified. 
    parser.add_argument('--organisms', default = "", nargs = "+") # list of organisms to check BUSCO completeness on
    parser.add_argument('--taxonomy_organisms', default = "", nargs = "+") # taxonomic level of organisms specified in organisms tag

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
    DATABASE_DIR = os.path.join(args.database)
    SAMPLE_DIR = args.sample_dir
    REF_FASTA = args.ref_fasta
    ## Convert reference FASTA file or directory to the proper path.
    if os.path.isdir(str(os.path.join(REF_FASTA))): # if reference fasta variable is a directory
        fasta_list = os.listdir(REF_FASTA)
        REFERENCE_FASTAS = [os.path.join(REFERENCE_DIR, curr) for curr in fasta_list]
        REF_FASTA = "combined_fastas"
    elif os.path.isfile(os.path.join(REFERENCE_DIR, REF_FASTA)):
        REFERENCE_FASTAS = [os.path.join(REFERENCE_DIR, REF_FASTA)]
    else:
        print("You need to either provide a single fasta reference file, or the name of a directory containing multiple reference FASTA files.")
        sys.exit(1)

    TAX_TAB = os.path.join(REFERENCE_DIR, args.tax_table)
    PROT_TAB = os.path.join(REFERENCE_DIR, args.protein_map)
    ALIGNMENT_CHOICE = args.alignment_choice
    IFPARALLEL = "series"
    if args.p:
        IFPARALLEL = "parallel"
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
    if args.individual_or_summary == "individual":
        if (BUSCO_FILE != "") & (os.path.isfile(BUSCO_FILE)):
            busco_file_read = read.csv(BUSCO_FILE, sep = "\t")
            ORGANISMS = list(busco_file_read.iloc[:,0])
            ORGANISMS_TAXONOMY = list(busco_file_read.iloc[:,1])
            print("Organisms and their taxonomy levels for BUSCO analysis were read from file.")
        else:
            print("No BUSCO file specified/found; using argument-specified organisms and taxonomy for BUSCO analysis.")

        if (len(ORGANISMS) != len(ORGANISMS_TAXONOMY)):
            print("Organisms and taxonomic specifications for BUSCO analysis do not contain the same number of entries. Please revise such that each organism flagged for BUSCO analysis also includes its original taxonomic level.")
            sys.exit(1)

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
    setupEukulele()

    ## Download the reference database if specified.
    if (args.reference_dir == "") | (not os.path.isdir(args.referencedir)) | \
        (not os.path.isdir(os.path.join(args.referencedir), REF_FASTA)):
        print("Automatically downloading database " + args.database.lower() + ". If you intended to " + \
              "use an existing database folder, be sure a reference FASTA, protein map, and taxonomy table " + \
              "are provided. Check the documentation for details.")
        REF_FASTA, TAX_TAB, PROT_TAB = downloadDatabase(args.database.lower(), args.alignment_choice)

    if SETUP:
        setupDatabases(create_tax_table, TAX_TAB, PROT_TAB, REF_FASTA, args.column, args.delimiter, original_tax_table, args.taxonomy_col_id, args.reformat, args.database, args.alignment_choice, REFERENCE_FASTAS, DATABASE_DIR, RERUN_RULES, args.strain_col_id, OUTPUTDIR)

    if (mets_or_mags == "mets"):
        ## Now for some TransDecoding ##
        samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == NT_EXT]
        print("Performing TransDecoder steps...", flush=True)
        if len(samples) == 0:
            print("No samples found in sample directory with specified nucleotide extension.")
        if ALIGNMENT:
            n_jobs_align = min(multiprocessing.cpu_count(), len(samples))
            transdecoder_res = Parallel(n_jobs=n_jobs_align)(delayed(transdecodeToPeptide)(sample_name) for sample_name in samples)
            all_codes = sum(transdecoder_res)
            if all_codes > 0:
                print("TransDecoder did not complete successfully; check log folder for details.")
                sys.exit(1)
            rcodes = [os.remove(curr) for curr in glob.glob("pipeliner*")]
    else:
        samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == PEP_EXT]

        if len(samples) == 0:
            print("No samples found in sample directory with specified peptide extension.")

    if ALIGNMENT:
        print("Performing alignment steps...", flush=True)
        ## First, we need to create a diamond database ##
        os.system("mkdir -p " + os.path.join(DATABASE_DIR, "diamond"))
        outfile_dmnd = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if args.alignment_choice == "blast":
            os.system("makeblastdb -in " + os.path.join(OUTPUTDIR, "concatfasta.fasta") + " -parse_seqids -blastdb_version " + str(5) + " -title " + str(database) + " -dbtype prot -out " + outfile_dmnd + " 1> log/db_create.err 2> log/db_create.log")
        else:
            os.system("diamond makedb --in " + os.path.join(OUTPUTDIR, "concatfasta.fasta") + " --db " + outfile_dmnd + " 1> log/db_create.err 2> log/db_create.log")
        
        ## Next to align against our database of choice ##
        n_jobs_align = min(multiprocessing.cpu_count(), len(samples))
        alignment_res = Parallel(n_jobs=n_jobs_align, prefer="threads")(delayed(align_to_database)(args.alignment_choice, sample_name, args.filter_metric, OUTPUTDIR, REF_FASTA, mets_or_mags, DATABASE_DIR, SAMPLE_DIR, RERUN_RULES, NT_EXT, PEP_EXT) for sample_name in samples)
        if any([((curr == None) | (curr == 1)) for curr in alignment_res]):
            print("Alignment did not complete successfully.")
            sys.exit(1)

        ## Next to do taxonomy estimation (NOTE: currently only supported with a config file). ##
        if (USE_SALMON_COUNTS == 1):
            if args.config_file != "":
                rc1 = namesToReads(args.config_file)

        print("Performing taxonomic estimation steps...", flush=True)
        outfiles = [os.path.join(OUTPUTDIR, mets_or_mags, samp + "-estimated-taxonomy.out") for samp in samples]
        n_jobs_align = min(multiprocessing.cpu_count(), len(alignment_res))
        for t in range(len(alignment_res)): 
            curr_out = place_taxonomy(TAX_TAB,args.cutoff_file,CONSENSUS_CUTOFF,\
                                                    PROT_TAB,USE_SALMON_COUNTS,NAMES_TO_READS,alignment_res[t],outfiles[t],\
                                                    IFPARALLEL,RERUN_RULES)

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
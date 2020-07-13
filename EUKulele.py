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
sys.path.insert(1, 'scripts')

import tax_placement
from tax_placement import *

#from query_busco import *

import visualize_results
from visualize_results import *

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

## FUNCTIONS TO USE IN PIPELINE ##

def setup():
    if args.create_tax_table:
        if args.original_tax_table == "":
            print("You must provide a taxonomy table via the argument 'original_tax_table' if you wish to run a taxonomy.")
            sys.exit(1)
        eukprot = ""
        if args.database == "eukprot":
            eukprot = " --euk-prot "
        rc1 = os.system("python scripts/create_protein_table --infile_peptide " + REF_FASTA + " --infile_taxonomy " + args.original_tax_table + " --output " + TAX_TAB  + " --outfile_json " + PROT_TAB + " --delim " + args.delimiter + " --strain_col_id " + args.strain_col_id + " --taxonomy_col_id " + args.taxonomy_col_id + " --column " + args.column + " --reformat_tax " + args.reformat + eukprot)

    ## Concatenate potential list of input FASTA files ##
    concatenated_file = os.path.join(OUTPUTDIR, "concatfasta.fasta")
    if (not os.path.isfile(concatenated_file)) | RERUN_RULES:
        space_delim = " ".join(REFERENCE_FASTAS)
        p1 = os.system("for currfile in " + space_delim + "; do ((cat $currfile | sed 's/\./N/g'); echo; echo) >> " + concatenated_file + "; done")
    else:
        print("Concatenated file already found in output directory; will not re-run step.", flush = True)

    output_log = "alignment_out.log"
    error_log = "alignment_err.log"
    if args.alignment_choice == "diamond":
        align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if (not os.path.isfile(align_db)) | RERUN_RULES:
            ## DIAMOND database creation ##
            db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa'))
            os.system("diamond makedb --in " + concatenated_file + " --db " + db + " 1> " + output_log + " 2> " + error_log)
        else:
            print("Diamond database file already created; will not re-create database.", flush = True)
    else:
        db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
        db_type = "prot"
        blast_version = 5
        os.system("makeblastdb -in " + concatenated_file + " -parse_seqids -blastdb_version " + str(blast_version) + " -title " + args.database + " -dbtype " + db_type + " -out " + db)

def transdecode_to_peptide(sample_name):
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder"))
    if (os.path.isfile(os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT))) | RERUN_RULES:
        print("TransDecoder file already detected for sample " + str(sample_name) + "; will not re-run step.")
        return 0
    returncode1 = os.system("TransDecoder.LongOrfs -t " + os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT) + " -m " + str(TRANSDECODERORFSIZE) + " 2> " + os.path.join("log", "transdecoder_error_" + sample_name + ".err") + " 1> " + os.path.join("log", "transdecoder_out_" + sample_name + ".out"))
    rc1 = returncode1
    returncode2 = os.system("TransDecoder.Predict -t " + os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT) + " --no_refine_starts 2>> " + os.path.join("log", "transdecoder_error_" + sample_name + ".err") + " 1>> " + os.path.join("log", "transdecoder_out_" + sample_name + ".out"))
    rc2 = returncode2
    if (rc1 + rc2) != 0:
        print("TransDecoder did not complete successfully for sample " + str(sample_name) + ". Check log/ folder for details.")
        sys.exit(1)
    merged_name = sample_name + "." + NT_EXT
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags))
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder"))
    os.replace(merged_name + ".transdecoder.pep", os.path.join(OUTPUTDIR, mets_or_mags,  sample_name + "." + PEP_EXT))
    os.replace(merged_name + ".transdecoder.cds", os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder", sample_name + ".fasta.transdecoder.cds"))
    os.replace(merged_name + ".transdecoder.gff3", os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder", sample_name + ".fasta.transdecoder.gff3"))
    os.replace(merged_name + ".transdecoder.bed", os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder", sample_name + ".fasta.transdecoder.bed"))
    shutil.rmtree(merged_name + ".transdecoder_dir*")
    return rc1 + rc2

def align_to_database(alignment_choice, sample_name, filter_metric):
    if alignment_choice == "diamond":
        os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "diamond"))
        diamond_out = os.path.join(OUTPUTDIR, mets_or_mags, "diamond", sample_name + ".diamond.out")
        if (os.path.isfile(diamond_out)) | RERUN_RULES:
            print("Diamond alignment file already detected; will not re-run step.")
            return diamond_out
        
        align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if mets_or_mags == "mets":
            fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
        else:
            fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)
        other = "--outfmt 6 -k 100 -e 1e-5"
        outfmt = 6
        k = 100
        e = 1e-5
        bitscore = 50
        diamond_log = os.path.join("log","diamond_align_" + sample_name + ".log")
        diamond_err = os.path.join("log","diamond_align_" + sample_name + ".err")
        if filter_metric == "bitscore":
            rc1 = os.system("diamond blastp --db " + align_db + " -q " + fasta + " -o " + diamond_out + " --outfmt " + str(outfmt) + " -k " + str(k) + " --min-score " + str(bitscore) + " 2> " + diamond_log + " 1> " + diamond_err)
        else:
            rc1 = os.system("diamond blastp --db " + align_db + " -q " + fasta + " -o " + diamond_out + " --outfmt " + str(outfmt) + " -k " + str(k) + " -e " + str(e) + " 2> " + diamond_log + " 1> " + diamond_err)
        if rc1 != 0:
            print("Diamond did not complete successfully.")
            os.system("rm -f " + diamond_out)
            return 1
        return diamond_out
    else:
        blast_out = os.path.join(OUTPUTDIR, mets_or_mags, "blast", sample_name + ".blast.txt")
        if (os.path.isfile(blast_out)) | RERUN_RULES:
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

def assign_taxonomy(sample_name):
    taxfile = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "-estimated-taxonomy.out")
    levels_directory = os.path.join(OUTPUTDIR, args.mets_or_mags, "levels"),
    max_dir = os.path.join(OUTPUTDIR, args.mets_or_mags)
    error_log = os.path.join("log", "tax_assign_" + sample_name + ".err")
    out_log = os.path.join("log", "tax_assign_" + sample_name + ".out")
    rc = os.system("python scripts/mag-stats.py --estimated-taxonomy-file " + taxfile + " --out-prefix " + sample_name + " --outdir " + levels_directory + " --max-out-dir " + max_dir + " 2> " + error_log + " 1> " + out_log)
    return rc

    
def run_busco(sample_name, outputdir, busco_db):
    if mets_or_mags == "mets":
        fastaname = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
    else:
        fastaname = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)
    os.system("chmod 755 scripts/run_busco.sh")
    rc1 = os.system("./scripts/run_busco.sh " + str(outputdir) + " static/busco_config.ini " + outputdir + "/config_" + sample_name + ".ini " + fastaname + " " + str(CPUS) + " " + busco_db)
    return rc1

parser = argparse.ArgumentParser()
parser.add_argument('subroutine', metavar="subroutine", nargs='?', type=str, default="all", choices = ["all","setup","alignment","busco"], help='Choice of subroutine to run.')
parser.add_argument('--mets_or_mags', required = True) 
parser.add_argument('--nucleotide_extension', default = ".fasta") 
parser.add_argument('--protein_extension', default = ".faa") 
parser.add_argument('--force_rerun', action='store_true', default=False)
parser.add_argument('--scratch', default = '../scratch') # the scratch location to store intermediate files

## SALMON OPTIONS ##
parser.add_argument('--use_salmon_counts', type = int, default = 0)
parser.add_argument('--salmon_dir') # salmon directory is required if use_salmon_counts is true.
parser.add_argument('--names_to_reads',default=0) # a file to be created or used if it exists that relates transcript names to salmon counts from the salmon directory 

## WHERE FILES ARE LOCATED ##
parser.add_argument('--database', default="mmetsp") # the name of the database to be used to assess the reads
parser.add_argument('--reference_dir', required = True) # folder containing the reference files for the chosen database
parser.add_argument('-o','--out_dir', dest = "out_dir", default = "output") # folder where the output will be written
parser.add_argument('--sample_dir', required = True) # folder where the input data is located (the protein or peptide files to be assessed)
# NOTE: make sure that the default reference.fasta name is used when later I provide a mechanism for downloading and formatting select databases.
parser.add_argument('--ref_fasta', default = "reference.pep.fa") # either a file in the reference directory where the fasta file for the database is located, or a directory containing multiple fasta files that constitute the database.
#parser.add_argument('--ref_fasta_ext', default = ".fasta") # if a directory is given for ref_fasta and the extension of the files differs from .fasta, specify it via this argument.

## TAXONOMY TABLE AND PROTEIN JSON FILE ## 
parser.add_argument('--create_tax_table', action='store_true') # include this file if you wish for the protein dictionary file and the taxonomy table to be created from the reference fasta(s) that you have provided. Otherwise, these files should be called "tax-table.txt" and "protein-map.json" and should be located in your reference_dir, unless specified via the below flags.
parser.add_argument('--original_tax_table', default = "", type = str) # this is required if you have specified "create_tax_table"
parser.add_argument('--strain_col_id',  type=str, default = 'strain_name') # the column which indicates the name of the strain in the taxonomy file
parser.add_argument('--taxonomy_col_id',  type=str, default = 'taxonomy') # the column which indicates the taxonomy of the strain in the taxonomy file
parser.add_argument('--column', type=str, default = 'SOURCE_ID') # can be numeric, zero-indexed, if it's a delimited part of header
# set to true if there is a column called "taxonomy" that we wish to split
parser.add_argument('--reformat_tax', dest='reformat', default=False, action='store_true') 

parser.add_argument('--delimiter', default = "/", type = str)
parser.add_argument('--tax_table', default = "tax-table.txt")
parser.add_argument('--protein_map', default = "protein-map.json")

## ALIGNMENT OPTIONS ##
parser.add_argument('--alignment_choice', default = "diamond", choices = ["diamond", "blast"])

## OPTIONS FOR CHECKING BUSCO COMPLETENESS FOR TAXONOMY ##
parser.add_argument('--busco_file', default = "", type = str) # if specified, the following two arguments ("--organisms" and "--taxonomy_organisms" are overwritten by the two columns of this tab-separated file
parser.add_argument('--individual_or_summary','-i',default="summary",choices=["summary","individual"])
# These arguments are used if "individual" is specified. 
parser.add_argument('--organisms', default = "", type = list, nargs = "+") # list of organisms to check BUSCO completeness on
parser.add_argument('--taxonomy_organisms', default = "", type = list, nargs = "+") # taxonomic level of organisms specified in organisms tag

## OTHER USER CHOICES ## 
parser.add_argument('--cutoff_file', default = "static/tax-cutoffs.yaml")
parser.add_argument('--filter_metric', default = "evalue", choices = ["pid", "evalue", "bitscore"])
parser.add_argument('--consensus_cutoff', default = 0.75, type = float)
parser.add_argument('--transdecoder_orfsize', default = 100, type = int)

parser.add_argument('--CPUs', default=1)
parser.add_argument('-p', action='store_true') # whether to run in parallel
parser.add_argument('--busco_threshold', default=50)

args = parser.parse_args()

## VARIABLES ##
CONSENSUS_CUTOFF = args.consensus_cutoff
REFERENCE_DIR = args.reference_dir
OUTPUTDIR = args.out_dir
DATABASE_DIR = os.path.join(REFERENCE_DIR, "database")
SAMPLE_DIR = args.sample_dir
REF_FASTA = args.ref_fasta
## Convert reference FASTA file or directory to the proper path.
if os.path.isdir(str(os.path.join(REF_FASTA))): # if reference fasta variable is a directory
    fasta_list = os.listdir(REF_FASTA)
    REFERENCE_FASTAS = [os.path.join(REFERENCE_DIR, curr) for curr in fasta_list]
    REF_FASTA = "combined_fastas"
    print(REF_FASTA)
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
mets_or_mags=args.mets_or_mags.lower()
if (mets_or_mags != "mets") & (mets_or_mags != "mags"):
    print("Only METs or MAGs are supported as input data types. Please update the 'mets_or_mags' flag accordingly.")
    sys.exit(1)
    
USE_SALMON_COUNTS = args.use_salmon_counts
SALMON_DIR = args.salmon_dir
NAMES_TO_READS = os.path.join(REFERENCE_DIR, args.names_to_reads)
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

## SETUP STEPS ##
print("Setting things up...")
#rc1 = os.system("conda activate EUKulele")
rc1 = subprocess.call(["activate", "euk-env/EUKulele"])
print(rc1)

if (rc1 != 0): # & (not os.path.isdir("./euk-env")):
    print("No EUKulele conda environment found; generating environment from EUKulele-env.yaml file...")
    os.system("conda env create -f EUKulele-env.yaml --force --prefix ./euk-env")
    p1 = subprocess.call(["activate", "EUKulele"])
    p1.wait()
    rc1 = subprocess.returncode
    if (rc != 0):
        print("Could not successfully generate and activate conda environment.")
        sys.exit(1)
        
os.system("mkdir -p " + OUTPUTDIR)
os.system("mkdir -p log")
if SETUP:
    setup()

if (mets_or_mags == "mets"):
    ## Now for some TransDecoding ##
    samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == NT_EXT]
    print("Performing TransDecoder steps...", flush=True)
    if len(samples) == 0:
        print("No samples found in sample directory with specified nucleotide extension.")
    if ALIGNMENT:
        n_jobs_align = min(multiprocessing.cpu_count(), len(samples))
        transdecoder_res = Parallel(n_jobs=n_jobs_align)(delayed(transdecode_to_peptide)(sample_name) for sample_name in samples)
        all_codes = sum(transdecoder_res)
        if all_codes > 0:
            print("TransDecoder did not complete successfully; check log folder for details.")
            sys.exit(1)
        rcodes = [os.remove(curr) for curr in glob.glob("pipeliner*")]
        #shutil.rmtree(merged_name + ".transdecoder_dir*")
else:
    samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == PEP_EXT]
    
    if len(samples) == 0:
        print("No samples found in sample directory with specified peptide extension.")
    
if ALIGNMENT:
    print("Performing alignment steps...", flush=True)
    ## Next to align against our database of choice ##
    n_jobs_align = min(multiprocessing.cpu_count(), len(samples))
    alignment_res = Parallel(n_jobs=n_jobs_align)(delayed(align_to_database)(args.alignment_choice, sample_name, args.filter_metric) for sample_name in samples)
    if any([(curr == None) for curr in alignment_res]):
        print("Alignment did not complete successfully.")
        sys.exit(1)

    ## Next to do taxonomy estimation ##
    if (USE_SALMON_COUNTS == 1):
        rc1 = os.system("python scripts/names_to_reads.py")

    print("Performing taxonomic estimation steps...", flush=True)
    outfiles = [os.path.join(OUTPUTDIR, mets_or_mags, samp + "-estimated-taxonomy.out") for samp in samples]
    n_jobs_align = min(multiprocessing.cpu_count(), len(alignment_res))
    taxonomy_res = Parallel(n_jobs=n_jobs_align)(delayed(place_taxonomy)(TAX_TAB,args.cutoff_file,CONSENSUS_CUTOFF,PROT_TAB,USE_SALMON_COUNTS,NAMES_TO_READS,alignment_res[t],outfiles[t],IFPARALLEL,RERUN_RULES) for t in range(len(alignment_res)))

    ## Now to visualize the taxonomy ##
    print("Performing taxonomic visualization steps...", flush=True)
    out_prefix = OUTPUTDIR.split("/")[-1]
    visualize_all_results(out_prefix, OUTPUTDIR, os.path.join(OUTPUTDIR, mets_or_mags), os.path.join(OUTPUTDIR, mets_or_mags), PEP_EXT, NT_EXT, USE_SALMON_COUNTS, RERUN_RULES)
    
    if mets_or_mags == "mags":
        print("Performing taxonomic assignment steps...", flush=True)
        n_jobs_viz = min(multiprocessing.cpu_count(), len(samples))
        assign_res = Parallel(n_jobs=n_jobs_viz)(delayed(assign_taxonomy)(samp) for samp in samples)
        if sum(assign_res) != 0:
            print("Taxonomic assignment of MAGs did not complete successfully. Check log files for details.")
            sys.exit(1)

print("Performing BUSCO steps...", flush=True)
if BUSCO:
    print("running busco")
    ## Run BUSCO on the full dataset ##
    busco_db = "eukaryota_odb10"
    busco_res = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(run_busco)(sample_name, os.path.join(OUTPUTDIR, "busco"), busco_db) for sample_name in samples)
    all_codes = sum(busco_res)
    if all_codes > 0:
        print("BUSCO did not complete successfully.")
        sys.exit(1)
    
## Assess BUSCO completeness on the most prevalent members of the metatranscriptome at each taxonomic level ##
if args.individual_or_summary == "individual":
    for sample_name in MTS:
        busco_table = os.path.join(OUTPUTDIR, "busco", sample_name, "full_table.tsv") # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
        taxtfile_stub = os.path.join(OUTPUTDIR,OUTPUTDIR.split("/")[-1]) # the prefix to specify where the taxonomy estimation output files are located
        
        if mets_or_mags == "mets":
            fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
        else:
            fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)
            
        query_busco_log = os.path.join("log","busco_query_" + sample_name + ".log")
        rc = os.system("python scripts/query_busco.py --organism_group " + " ".join(args.organisms) + " --taxonomic_level " + " ".join(args.taxonomy_organisms) + " --output_dir " + OUTPUTDIR + " --fasta_file " + fasta + " --sample_name " + sample_name + " --taxonomy_file_prefix " + taxfile_stub + " --tax_table " + TAX_TAB + " --busco_out " + busco_table + "-i individual > " + query_busco_log)
        if rc != 0:
            print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")
else:
    for sample_name in MTS:
        busco_table = os.path.join(OUTPUTDIR, "busco", sample_name, "full_table.tsv") # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
        taxtfile_stub = os.path.join(OUTPUTDIR,OUTPUTDIR.split("/")[-1]) # the prefix to specify where the taxonomy estimation output files are located
        
        if mets_or_mags == "mets":
            fasta = os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT) 
        else:
            fasta = os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT)
            
        query_busco_log = os.path.join("log","busco_query_" + sample_name + ".log")
        rc = os.system("python scripts/query_busco.py --output_dir " + OUTPUTDIR + " --fasta_file " + fasta + " --sample_name " + sample_name + " --taxonomy_file_prefix " + taxfile_stub + " --tax_table " + TAX_TAB + " --busco_out " + busco_table + " -i summary > " + query_busco_log)
        if rc != 0:
            print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")

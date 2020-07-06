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
from joblib import Parallel, delayed
sys.path.insert(1, 'scripts')

import tax_placement
from tax_placement import *

import query_busco
from query_busco import *

import visualize_results
from visualize_results import *

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

parser = argparse.ArgumentParser()
parser.add_argument('--mets_or_mags', required = True) 
parser.add_argument('--nucleotide_extension', default = ".fasta") 
parser.add_argument('--protein_extension', default = ".faa") 
parser.add_argument('--scratch', default = '../scratch') # the scratch location to store intermediate files

## SALMON OPTIONS ##
parser.add_argument('--use_salmon_counts', type = int, default = 0)
parser.add_argument('--salmon_dir') # salmon directory is required if use_salmon_counts is true.
parser.add_argument('--names_to_reads') # a file to be created or used if it exists that relates transcript names to salmon counts from the salmon directory 

## WHERE FILES ARE LOCATED ##
parser.add_argument('--database', required = True) # the name of the database to be used to assess the reads
parser.add_argument('--reference_dir', required = True) # folder containing the reference files for the chosen database
parser.add_argument('-o','--out_dir', dest = "out_dir", required = True) # folder where the output will be written
parser.add_argument('--sample_dir', required = True) # folder where the input data is located (the protein or peptide files to be assessed)
parser.add_argument('--ref_fasta', required = True) # either a file in the reference directory where the fasta file for the database is located, or a directory containing multiple fasta files that constitute the database.
parser.add_argument('--ref_fasta_ext', default = ".fasta") # if a directory is given for ref_fasta and the extension of the files differs from .fasta, specify it via this argument.

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

## OTHER USER CHOICES ## 
parser.add_argument('--cutoff', default = "static/tax-cutoffs.yaml")
parser.add_argument('--cutoff_metric', default = "pid", choices = ["pid", "evalue"])
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
TAX_TAB = os.path.join(REFERENCE_DIR, args.tax_table)
PROT_TAB = os.path.join(REFERENCE_DIR, args.protein_map)
ALIGNMENT_CHOICE = config['alignment_choice']
IFPARALLEL = config['choose_parallel']
OUTPUT_EXTENSION = "txt"
DBEXTENSION = ""
TRANSDECODERORFSIZE=config['transdecoder_orfsize']
if ALIGNMENT_CHOICE == "diamond":
    OUTPUT_EXTENSION = "out"
    DBEXTENSION = ".dmnd"
NT_EXT = config['nucleotide_extension'].strip('.')
PEP_EXT = config['protein_extension'].strip('.')
mets_or_mags=config['mets_or_mags'].lower()
USE_SALMON_COUNTS = config["use_salmon_counts"]
SALMON_DIR = config["salmon_dir"]
NAMES_TO_READS = config["names_to_reads"]

## SETUP STEPS ##
if args.create_tax_table:
    if args.original_tax_table == "":
        print("You must provide a taxonomy table via the argument 'original_tax_table' if you wish to run a taxonomy.")
        sys.exit(1)
    eukprot = ""
    if args.database == "eukprot":
        eukprot = " --euk-prot "
    p1 = subprocess.Popen(["python scripts/create_protein_table --infile_peptide " + REF_FASTA + " --infile_taxonomy " + args.original_tax_table + " --output " + TAX_TAB  + " --outfile_json " + PROT_TAB + " --delim " + args.delimiter + " --strain_col_id " + args.strain_col_id + " --taxonomy_col_id " + args.taxonomy_col_id + " --column " + args.column + " --reformat_tax " + args.reformat + eukprot])
    p1.wait()
    
## Concatenate potential list of input FASTA files ##
space_delim = " ".join(REFERENCE_FASTAS)
concatenated_file = os.path.join(OUTPUTDIR, "concatfasta.fasta")
p1 = subprocess.Popen(["for currfile in " + space_delim + "; do ((cat $currfile | sed 's/\./N/g'); echo; echo) >> " + concatenated_file + "; done"])
p1.wait()

if args.alignment_choice == "diamond":
    ## DIAMOND database creation ##
    db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa'))
    p1 = subprocess.Popen(["diamond makedb --in " + concatenated_file + " --db " + db])
    p1.wait()
else:
    db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
    db_type = "prot"
    blast_version = 5
    p1 = subprocess.Popen(["makeblastdb -in " + concatenated_file + " -parse_seqids -blastdb_version " + str(blast_version) + " -title ". + args.database + " -dbtype " + db_type + " -out " + db])
    p1.wait()

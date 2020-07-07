#!/usr/bin/env python

# Query BUSCO for a given species or functional group. 
# A script that takes a given functional group as input, along with the taxonomic level of that group
# (e.g. Phaeocystis antarctica, species), and then checks for BUSCO completeness among the contigs
# identified as that taxonomic level or lower, also evaluating the number of copies of the BUSCO
# matches to differentiate between multiple strains.

import pandas as pd
import numpy as np
import os
import sys

__author__ = "Arianna Krinos, Harriet Alexander"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

parser = argparse.ArgumentParser()
parser.add_argument('--organism_group') # the focal name of species/genus/order/class etc.
parser.add_argument('--taxonomic_level') # the taxonomic level of the specified focal name.
parser.add_argument('--fasta_file') # the fasta file from which we pull sequences for the mock transcriptome.
parser.add_argument('--taxonomy_file_prefix') # the taxonomy file we use to create the mock transcriptome.
parser.add_argument('--tax_table') # the taxonomy table to get the full classification of the organism as assessed by the database being used.
parser.add_argument('--sample_name') # name of the original sample being assessed.
parser.add_argument('--download_busco',action='store_true') # if specified, we download BUSCO file from the url in the next argument
parser.add_argument('--busco_url',default=0)
parser.add_argument('--busco_location',default="busco") # location to store the BUSCO tar reference
parser.add_argument('--output_dir',default="output")
parser.add_argument('--available_cpus',default=1)
parser.add_argument('--busco_threshold',default=50)

def read_in_taxonomy(infile):
    with open(infile, 'rb') as f:
        result = chardet.detect(f.read())
    tax_out = pd.read_csv(infile, sep='\t',encoding=result['encoding'])
    tax_out.columns = tax_out.columns.str.lower()
    tax_out = tax_out.set_index('source_id')
    return tax_out

args = parser.parse_args()
organism = args.organism_group
taxonomy = args.taxonomic_level
tax_table = read_in_taxonomy(args.tax_table)
full_taxonomy = tax_table.loc[tax_table[taxonomy] == organism] #[organism in curr for curr in tax_table[taxonomy]]
if len(full_taxonomy).index < 1:
    print("No taxonomy found for that organism and taxonomic level.")
    sys.exit(1)
level_hierarchy = ['supergroup','division','class','order','family','genus','species']
curr_level = [ind for ind in range(len(level_hierarchy)) if level_hierarchy[ind] == taxonomy]
max_level = len(level_hierarchy) - 1

success = 0
busco_scores = []
while (curr_level >= 0):
    
    #### GET THE CURRENT LEVEL OF TAXONOMY FROM THE TAX TABLE FILE ####
    curr_taxonomy = str(full_taxonomy[level_hierarchy[curr_level]][0])
    if (curr_taxonomy == "") | (curr_taxonomy.lower() == "nan"):
        print("No taxonomy found at level " + level_hierarchy[curr_level])
        continue
        
    #### CREATE THE MOCK TRANSCRIPTOME BY PULLING BY TAXONOMIC LEVEL ####
    taxonomy_file = pd.read_csv(args.taxonomy_file_prefix + "_all_" + str(level_hierarchy[curr_level]) + "_counts.csv", sep="\t")
    taxonomy_file = taxonomy_file.loc[taxonomy_file[level_hierarchy[curr_level].capitalize()] == curr_taxonomy]
    transcripts_to_search = list(taxonomy_file["GroupedTranscripts"])
    
    ## Write the candidate transcripts to a file for easy grepping ##
    searchfile = organism + "_" + args.taxonomic_level + "_" + args.sample_name + "_transcriptnames.txt"
    with open(searchfile, 'w') as filehandle:
        for transcript_name in transcripts_to_search:
            filehandle.write(transcript_name + '\n')
    
    mock_file_name = organism + "_" + args.taxonomic_level + "_" + args.sample_name + "_curr.fasta"
    os.system("grep -w -A 2 -f " + searchfile + " " + args.fasta_file + " --no-group-separator > " + mock_file_name)

    #### USING BUSCO TO ASSESS COMPLETENESS OF THE MOCK TRANSCRIPTOME ####
    mock_file = mock_file_name # here will go path to FASTA file we create and assess completeness against

    if (args.download_busco) & (busco_url != 0):
        os.system("wget -O busco_db.tar.gz " + args.busco_url)
        os.system("mkdir -p " + args.busco_location)
        os.system("tar -xzf busco_db.tar.gz -C " + args.busco_location)
    elif (args.download_busco):
        print("You asked to download a BUSCO database, but didn't provide a URL for one.")
        sys.exit(1)

    ## By default, BUSCO output will just be stored below the output directory
    busco_loc = os.path.join(args.output_dir, "busco_run_" + organism + "_" + taxonomy)
    busco_db_name = "eukaryota_odb10" # we can also change this to our downloaded BUSCO file; assess this in the future
    os.system("cp " + os.path.join("..","static","busco_config.ini") + " " + os.path.join(busco_loc,"config.ini"))
    os.system("sed -i '/out = /c\out = " + organism + "' " + os.path.join(busco_loc,"config.ini")) # the name of the output files
    os.system("sed -i '/out_path = /c\out_path = " + busco_loc + "' " + os.path.join(busco_loc,"config.ini")) # what directory the output will be stored in
    os.system("busco -i " + mock_file + " -l " + busco_db_name + " -m transcriptome --cpu " + str(args.available_cpus) + " --config " + os.path.join(busco_loc,"config.ini") + " -f")

    #### PROCESS BUSCO OUTPUT ####
    busco_short_result = glob.glob(os.path.join(busco_loc,"short_summary*.txt"))
    print(busco_short_result)
    with open(busco_short_result[0], 'r') as file:
        data = file.read().replace('\n', '')
    busco_completeness = float(data.split("C:")[1].split("%")[0])
    busco_scores.append(busco_completeness)
    if busco_completeness >= args.busco_threshold:
        success = 1
        break
    curr_level = curr_level - 1

report_file = os.path.join(OUTPUTDIR, "busco_run", organism, args.taxonomic_level, args.sample_name + "_report.txt")
#report_file = os.path.join(args.output_dir, "busco_run_" + organism + "_" + args.taxonomic_level + "_" + args.sample_name + "_report.txt")
reported = open(report_file,"w")
if success == 1:
    file_written = os.path.join(args.output_dir, organism + level_hierarchy[curr_level] + ".fasta")
    os.system("mv " + mock_file + " " + file_written)
    reported.write("Taxonomy file successfully completed with BUSCO completeness " + str(busco_completeness) + "% at location " + str(file_written) + ". \n")
    reported.write("The BUSCO scores found at the various taxonomic levels were: " + str(busco_scores))
else:
    reported.write("Sufficient BUSCO completeness not found at threshold " + str(args.busco_threshold) + "%. \n")
    reported.write("The BUSCO scores found at the various taxonomic levels were: " + str(busco_scores))
reported.close()
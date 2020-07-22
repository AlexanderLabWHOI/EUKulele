# Query BUSCO for a given species or functional group. 
# A script that takes a given functional group as input, along with the taxonomic level of that group
# (e.g. Phaeocystis antarctica, species), and then checks for BUSCO completeness among the contigs
# identified as that taxonomic level or lower, also evaluating the number of copies of the BUSCO
# matches to differentiate between multiple strains. # !/ usr/bin/env python3

# os.path.join("EUKulele/tests/aux_data/test_out","test_out")
# python query_busco.py --organism_group Chromera --taxonomic_level genus --output_dir check_dir --fasta_file  --sample_name sample_0 --taxonomy_file_prefix   --tax_table " + TAX_TAB + " --busco_out " + busco_table

import pandas as pd
import numpy as np
import os
import sys
import argparse
import chardet
import glob
import multiprocessing
from joblib import Parallel, delayed

__author__ = "Arianna Krinos, Harriet Alexander"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

parser = argparse.ArgumentParser()
parser.add_argument('--busco_out',default="busco") # the output from the BUSCO run on the full sample file
parser.add_argument('--individual_or_summary','-i',default="summary",choices=["summary","individual"])
# Not necessary if we are running in summary mode.
parser.add_argument('--organism_group', nargs = "+", default = []) # the focal name(s) of species/genus/order/class etc.
parser.add_argument('--taxonomic_level', nargs = "+", default = []) # the taxonomic level(s) of the specified focal name.
parser.add_argument('--fasta_file') # the fasta file from which we pull sequences for the mock transcriptome.
parser.add_argument('--taxonomy_file_prefix') # the taxonomy file we use to create the mock transcriptome.
parser.add_argument('--tax_table') # the taxonomy table to get the full classification of the organism as assessed by the database being used.
parser.add_argument('--sample_name') # name of the original sample being assessed.
parser.add_argument('--download_busco',action='store_true') # if specified, we download BUSCO file from the url in the next argument
parser.add_argument('--create_fasta',action='store_true') # if specified, we create a "transcriptome fasta" when we query for the organisms
parser.add_argument('--busco_url',default=0)
parser.add_argument('--busco_location',default="busco") # location to store the BUSCO tar reference
parser.add_argument('--output_dir',default="output")
parser.add_argument('--available_cpus',default=1)
parser.add_argument('--busco_threshold',default=50)
parser.add_argument('--write_transcript_file',default=False,action='store_true') # whether to write an actual file with the subsetted transcriptome

level_hierarchy = ['supergroup','division','class','order','family','genus','species']

def evaluate_organism(organism, taxonomy, tax_table, create_fasta, write_transcript_file):
    organism_format = " ".join(organism.split("_"))
    full_taxonomy = tax_table.loc[[(organism_format in curr) for curr in tax_table[taxonomy]],:]
    if len(full_taxonomy.index) < 1:
        print("No taxonomy found for that organism and taxonomic level.")
        sys.exit(1)
    curr_level = [ind for ind in range(len(level_hierarchy)) if level_hierarchy[ind] == taxonomy][0]
    max_level = len(level_hierarchy) - 1

    success = 0
    success_level = ""
    busco_scores = []
    levels_out = []
    percent_multiples = []
    number_duplicated = [] 
    number_tripled = []
    number_quadrupled = []
    number_higher_mult = []
    number_covered = []
    busco_out_file = pd.read_csv(args.busco_out, sep = "\t", comment = "#", names = ["BuscoID","Status","Sequence","Score","Length"])
    select_inds = [ (busco_out_file.Status[curr] == "Complete") | (busco_out_file.Status[curr] == "Fragmented") | (busco_out_file.Status[curr] == "Duplicated") for curr in range(len(busco_out_file.index))]
    good_buscos = busco_out_file.loc[select_inds,:]
    good_busco_sequences = [curr.split(".")[0] for curr in list(good_buscos.Sequence)]
    Covered_IDs = list(good_buscos.BuscoID)
    total_buscos = len(set(list(busco_out_file.BuscoID)))
    
    while (curr_level >= 0):

        #### GET THE CURRENT LEVEL OF TAXONOMY FROM THE TAX TABLE FILE ####
        curr_tax_list = set(list(full_taxonomy[level_hierarchy[curr_level]]))
        if len(curr_tax_list) > 1:
            print(curr_tax_list)
            print("More than 1 unique match found; using both matches: " + str("".join(curr_tax_list)))
                  #using the first match (" + str(list(curr_tax_list)[0]) + ")")
        curr_taxonomy = ";".join(curr_tax_list) #str(list(curr_tax_list)[0])
        if (curr_taxonomy == "") | (curr_taxonomy.lower() == "nan"):
            print("No taxonomy found at level " + level_hierarchy[curr_level])
            continue

        #### CREATE A "MOCK TRANSCRIPTOME" BY PULLING BY TAXONOMIC LEVEL ####
        taxonomy_file = pd.read_csv(args.taxonomy_file_prefix + "_all_" + str(level_hierarchy[curr_level]) + "_counts.csv", sep=",",header=0)
        taxonomy_file = taxonomy_file.loc[[tax in curr_taxonomy for tax in list(taxonomy_file[level_hierarchy[curr_level].capitalize()])],:]
        transcripts_to_search = list(taxonomy_file["GroupedTranscripts"])
        transcripts_to_search_sep = []
        for transcript_name in transcripts_to_search:
            transcripts_to_search_sep.extend([curr.split(".")[0] for curr in transcript_name.split(";")])

        set_transcripts_to_search = set(transcripts_to_search_sep)
        good_busco_sequences_list = list(good_busco_sequences)
        BUSCOs_covered_all = [Covered_IDs[curr] for curr in range(len(good_busco_sequences_list)) if good_busco_sequences_list[curr] in list(set_transcripts_to_search)]
        BUSCOs_covered = set(BUSCOs_covered_all)
        number_appearences = [BUSCOs_covered_all.count(curr_busco) for curr_busco in list(BUSCOs_covered)]
        multiples = [curr_appear for curr_appear in number_appearences if curr_appear >= 2]
        
        ## KEEP TRACK OF WHETHER DUPLICATES/TRIPLES/ETC. ARE COMMON ##
        number_covered.append(number_appearences)
        number_duplicated.append(number_appearences.count(2))
        number_tripled.append(number_appearences.count(3))
        number_quadrupled.append(number_appearences.count(4))
        number_higher_mult.append(len([curr_appear for curr_appear in number_appearences if curr_appear > 4]))
        prop_duplicated = 0
        if len(BUSCOs_covered) > 0:
            prop_duplicated = len(multiples) / len(BUSCOs_covered) * 100
            
        percent_multiples.append(prop_duplicated)
        
        busco_completeness = len(BUSCOs_covered) / total_buscos * 100
        busco_scores.append(busco_completeness)
        levels_out.append(level_hierarchy[curr_level])
        if busco_completeness >= args.busco_threshold:
            success = 1
            success_level = level_hierarchy[curr_level]
        curr_level = curr_level - 1

    report_dir = os.path.join(args.output_dir, organism, taxonomy)
    os.system("mkdir -p " + os.path.join(args.output_dir, organism))
    os.system("mkdir -p " + report_dir)

    report_file = os.path.join(report_dir, args.sample_name + "_report.txt")
    reported = open(report_file,"w")
    if success == 1:
        file_written = os.path.join(args.output_dir, organism, level_hierarchy[curr_level] + "_" + args.sample_name + ".txt")
        if (write_transcript_file):
            with open(file_written, 'w') as filehandle:
                for transcript_name in transcripts_to_search_sep:
                    filehandle.write(transcript_name + '\n')
        if (create_fasta):
            mock_file_name = organism + "_" + level_hierarchy[curr_level] + "_" + args.sample_name + "_BUSCO_complete.fasta"
            os.system("grep -w -A 2 -f " + file_written + " " + args.fasta_file + " --no-group-separator > " + mock_file_name)

        reported.write("Taxonomy file successfully completed with BUSCO completeness " + str(busco_completeness) + "% at location " + str(file_written) + "This was at taxonomic level " + str(success_level) + ". The file containing the transcript names for the mock transcriptome corresponding to this taxonomic level is located here: " + str(file_written) + ".\n")
        reported.write("The BUSCO scores found at the various taxonomic levels (Supergroup to " + str(args.taxonomic_level) + ") were: " + str(busco_scores))
    else:
        reported.write("Sufficient BUSCO completeness not found at threshold " + str(args.busco_threshold) + "%. \n")
        reported.write("The BUSCO scores found at the various taxonomic levels (Supergroup to " + str(args.taxonomic_level) + ") were: " + str(busco_scores.reverse()) + "\n")
    reported.close()
    return pd.DataFrame({"Organism":[organism] * len(levels_out),"TaxonomicLevel":levels_out,"BuscoCompleteness":busco_scores,"NumberCovered":number_covered,"CtTwoCopies":number_duplicated,"CtThreeCopies":number_tripled,"CtFourCopies":number_quadrupled,"CtFivePlusCopies":number_higher_mult,"PercentageDuplicated":percent_multiples})

def read_in_taxonomy(infile):
    with open(infile, 'rb') as f:
        result = chardet.detect(f.read())
    tax_out = pd.read_csv(infile, sep='\t',encoding=result['encoding'])
    tax_out.columns = tax_out.columns.str.lower()
    tax_out = tax_out.set_index('source_id')
    return tax_out

args = parser.parse_args()
#organism_format = (" ".join(args.organism_group) #" ".join(args.organism_group.split("_"))
organism = args.organism_group
taxonomy = args.taxonomic_level
tax_table = read_in_taxonomy(args.tax_table)
            
#full_taxonomy = tax_table.loc[[(organism_format in curr) for curr in tax_table[taxonomy]],:]
            
if (args.individual_or_summary == "individual") & ((len(args.organism_group) == 0) | (len(args.taxonomic_level) == 0)):
    print("You specified individual mode, but then did not provide a taxonomic group and/or accompanying taxonomic level.")
    sys.exit(1)
if (len(args.organism_group) == 0) != (len(args.taxonomic_level) == 0):
    print("The number of organisms you specified is not equal to the number of taxonomic levels you specified. " + str(len(args.organism_group)) + " organisms were specified, while " + str(args.taxonomic_level) + " taxonomic levels were specified.")
    sys.exit(1)

if (args.individual_or_summary == "individual"):
    results_frame = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(evaluate_organism)(organism[curr], taxonomy[curr], tax_table, args.create_fasta, args.write_transcript_file) for curr in range(len(organism)))
    results_frame = pd.concat(results_frame)
    os.system("mkdir -p " + os.path.join(args.output_dir, "busco_assessment", args.sample_name, "individual"))
    results_frame.to_csv(path_or_buf = os.path.join(args.output_dir, "busco_assessment", args.sample_name, "individual", "summary_" + args.sample_name + ".tsv"), sep = "\t")
else:
    for curr_level in level_hierarchy:
        taxonomy = level_hierarchy[curr_level]
        taxonomy_file = pd.read_csv(args.taxonomy_file_prefix + "_all_" + str(level_hierarchy[curr_level]) + "_counts.csv", sep=",",header=0)
        curr_frame = taxonomy_file.nlargest(multiprocessing.cpu_count(), 'NumTranscripts') # change this to number of CPUs 
        organisms = list(curr_frame[level_hierarchy[curr_level].capitalize()])
        # now do multiprocessing for each of these organisms and get BUSCO results for each
        # need to pull in argument for whether or not to save CSVs 
        results_frame = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(
    evaluate_organism)(organism, taxonomy, tax_table, args.create_fasta, args.write_transcript_file) for organism in range(len(organisms)))
        results_frame = pd.concat(results_frame)
        os.system("mkdir -p " + os.path.join(args.output_dir, "busco_assessment", args.sample_name, taxonomy + "_combined"))
        results_frame.to_csv(path_or_buf = os.path.join(args.output_dir, "busco_assessment", args.sample_name, taxonomy + "_combined", "summary_" + taxonomy + "_" +  args.sample_name + ".tsv"), sep = "\t")
    # here's where we're doing the cascade through the most popular. Requires that we create a helper function to grab then evaluate those.
    # TODO: assessing the copies of each of the BUSCOs found in the cascade. 
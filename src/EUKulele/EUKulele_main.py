'''
Software for taxonomic identification of eukaryotes.
'''

import os
import sys
import argparse
import multiprocessing
import pandas as pd
import datetime

from EUKulele.download_database import downloadDatabase
from EUKulele.manage_steps import manageEukulele
from EUKulele.busco_runner import readBuscoFile
from EUKulele.busco_runner import configRunBusco
from EUKulele.busco_runner import manageBuscoQuery

from scripts.names_to_reads import namesToReads

__author__ = "Harriet Alexander, Arianna Krinos"
__copyright__ = "EUKulele"
__license__ = "MIT"
__email__ = "akrinos@mit.edu"

def main(args_in):
    '''
    Main function for calling subfunctions and running EUKulele.
    '''

    parser = argparse.ArgumentParser(
        description='Thanks for using EUKulele! EUKulele is a standalone taxonomic '+\
                    'annotation software.\n EUKulele is designed primarily for marine '+\
                    'microbial eukaryotes. Check the README for further information.',
        usage='eukulele [subroutine] --mets_or_mags [dataset_type] --sample_dir ' +\
              '[sample_directory] --reference_dir [reference_database_location] ' +\
              '[all other options]')

    parser.add_argument('subroutine', metavar="subroutine", nargs='?', type=str,
                        default="all",
                        choices = ["","all","download","setup","alignment",
                                   "busco","coregenes"],
                        help='Choice of subroutine to run.')

    parser.add_argument('-v', '--version', dest = "version", default=False,
                        action='store_true')
    parser.add_argument('-m', '--mets_or_mags', dest = "mets_or_mags", required = False,
                        default = "")
    parser.add_argument('--n_ext', '--nucleotide_extension', dest = "nucleotide_extension",
                        default = ".fasta")
    parser.add_argument('--p_ext', '--protein_extension',
                        dest = "protein_extension",
                        default = ".faa")
    parser.add_argument('-f', '--force_rerun', action='store_true', default=False)
    parser.add_argument('--scratch', default = '../scratch',
                        help = "The scratch location to store intermediate files.")
    parser.add_argument('--config_file', default = '')
    parser.add_argument('--perc_mem', dest = "perc_mem",
                        default = 0.75,
                        help = "The percentage of the total available memory which should "+\
                               "be targeted for use by processes.")

    ## SALMON OPTIONS ##
    parser.add_argument('--use_salmon_counts', action='store_true', default=False)
    parser.add_argument('--salmon_dir',
                        help = "Salmon directory is required if use_salmon_counts is true.")
    parser.add_argument('--names_to_reads',default=0, help = "A file to be created or " +\
                        "used if it exists that relates transcript names to salmon counts "+\
                        "from the salmon directory.")

    ## WHERE FILES ARE LOCATED ##
    parser.add_argument('-d', '--database', default="marmmetsp",
                        help = "The name of the database to be used to assess the reads.")
    parser.add_argument('-o', '--out_dir', '--output_dir', dest = "out_dir", default = "output",
                        help = "Folder where the output will be written.")
    parser.add_argument('-s', '--sample_dir', required = False, dest = "sample_dir",
                        default = "nan",
                        help = "Folder where the input data is located (the protein or "+\
                               "peptide files to be assessed).")

    ## ONLY SPECIFY THESE ARGUMENTS IF YOU HAVE ALREADY PROVIDED AND FORMATTED YOUR OWN DATABASE ##
    parser.add_argument('--reference_dir', default=".",
                        help = "Folder containing the reference files for the chosen database.")
    parser.add_argument('--ref_fasta', default = "reference.pep.fa",
                        help = "Either a file in the reference directory where the " +\
                               "fasta file for the database is located, or a directory " +\
                               "containing multiple fasta files that " +\
                               "constitute the database.")

    parser.add_argument('--tax_table', default = "tax-table.txt")
    parser.add_argument('--protein_map', default = "prot-map.json")

    ## ALIGNMENT OPTIONS ##
    parser.add_argument('--alignment_choice', default = "diamond", choices = ["diamond", "blast"])

    ## OPTIONS FOR CHECKING BUSCO COMPLETENESS FOR TAXONOMY ##
    parser.add_argument('--busco_file', default = "", type = str,
                        help = "If specified, the following two arguments ('--organisms' and " +\
                        "'--taxonomy_organisms' are overwritten by the two columns of this "+\
                        "tab-separated file.")
    parser.add_argument('--individual_or_summary', default="summary",
                        choices=["summary","individual"],
                        help = "These arguments are used if 'individual' is specified.")
    parser.add_argument('-i', '--individual', dest = "individual_tag", action='store_true',
                        default=False)
    parser.add_argument('--organisms', default = "", nargs = "+",
                        help = "List of organisms to check BUSCO completeness on.")
    parser.add_argument('--taxonomy_organisms', default = "", nargs = "+",
                        help = "Taxonomic level of organisms specified in organisms tag.")

    ## OTHER USER CHOICES ##
    cutoff_file = "tax-cutoffs.yaml"
    parser.add_argument('--cutoff_file', default = cutoff_file)
    parser.add_argument('--filter_metric', default = "evalue",
                        choices = ["pid", "evalue", "bitscore"])
    parser.add_argument('--consensus_cutoff', default = 0.75, type = float)
    parser.add_argument('--transdecoder_orfsize', default = 100, type = int)

    parser.add_argument('--CPUs', default=multiprocessing.cpu_count())
    parser.add_argument('--busco_threshold', default=50)
    parser.add_argument('--create_fasta', action='store_true', default=False,
                       help = "Whether to create FASTA files containing ID'd transcripts "+\
                              "during BUSCO analysis.")
    parser.add_argument('--run_transdecoder', action='store_true', default=False,
                       help = "Whether TransDecoder should be run on metatranscriptomic samples. "+\
                        "Otherwise, BLASTp is run if protein translated samples are provided" +\
                        "otherwise BLASTx is run on nucleotide samples.")


    parser.add_argument('--test', action='store_true', default=False,
                       help = "Whether we're just running a test and should not execute downloads.")
  
    args = parser.parse_args(list(filter(None, args_in.split(" "))))
    if (args.mets_or_mags == "") & (args.subroutine != "download") & (not args.version):
        print("METs or MAGs argument (-m/--mets_or_mags) is required with one of 'mets' or 'mags'.")
        sys.exit(1)
    if (args.sample_dir == "nan") & (args.subroutine != "download") & (not args.version):
        print("A sample directory must be specified (-s/--sample_dir).")
        sys.exit(1)

    ## VARIABLES ##
    test_var = args.test
    consensus_cutoff = args.consensus_cutoff
    reference_dir = args.reference_dir
    output_dir = args.out_dir
    sample_dir = args.sample_dir
    ref_fasta = args.ref_fasta
    perc_mem = args.perc_mem

    alignment_choice = args.alignment_choice
    transdecoder_orf_size=args.transdecoder_orfsize
    nt_ext = args.nucleotide_extension.strip('.')
    pep_ext = args.protein_extension.strip('.')
    mets_or_mags = args.mets_or_mags.lower()
    individual_or_summary = args.individual_or_summary
    if args.individual_tag:
        individual_or_summary = "individual"

    if (mets_or_mags != "mets") & (mets_or_mags != "mags") & \
       (args.subroutine != "download") & (not args.version):
        print("Only METs or MAGs are supported as input data types. "+\
              "Please update the 'mets_or_mags' flag accordingly.")
        sys.exit(1)

    use_salmon_counts = args.use_salmon_counts
    salmon_dir = args.salmon_dir
    names_to_reads = os.path.join(reference_dir, str(args.names_to_reads))     

    organisms = args.organisms
    organisms_taxonomy = args.taxonomy_organisms
    busco_file = args.busco_file
    rerun_rules = args.force_rerun
    run_transdecoder = args.run_transdecoder

    organisms, organisms_taxonomy = readBuscoFile(individual_or_summary, busco_file,
                                                  organisms, organisms_taxonomy)

    download_choice = False
    setup_choice = False
    alignment_choice_select = False
    busco_choice = False
    core_genes = False

    if args.version:
        test_var = True
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "static", "VERSION")
        file_read = open(filename, "r")
        print("The current EUKulele version is",file_read.read())

    if args.subroutine == "download":
        download_choice = True
    if (args.subroutine == "all") | (args.subroutine == "setup"):
        setup_choice = True
    if (args.subroutine == "all") | (args.subroutine == "alignment"):
        alignment_choice_select = True
    if (args.subroutine == "all") | (args.subroutine == "busco"):
        busco_choice = True
    if (args.subroutine == "all") | (args.subroutine == "coregenes"):
        core_genes = True

    if test_var:
        setup_choice = False
        alignment_choice_select = False
        busco_choice = False
        core_genes = False
    else:
        ## SETUP STEPS / DOWNLOAD DEPENDENCIES ##
        if not download_choice:
            manageEukulele(piece = "setup_eukulele", output_dir = output_dir)
            samples = manageEukulele(piece = "get_samples", mets_or_mags = mets_or_mags,
                                     sample_dir = sample_dir, nt_ext = nt_ext,
                                     pep_ext = pep_ext)
        else:
            os.system("mkdir -p " + output_dir)
            os.system("mkdir -p " + os.path.join(output_dir, "log"))

        tax_tab = os.path.join(reference_dir, args.tax_table)
        prot_tab = os.path.join(reference_dir, args.protein_map)

        ## Download the reference database if specified. ##
        ## First, check if files are in originally specified reference directory ##
        if (not os.path.isfile(os.path.join(reference_dir, ref_fasta))) | \
           (not os.path.isfile(tax_tab)) | \
           (not os.path.isfile(prot_tab)):
            reference_dir = os.path.join(reference_dir, args.database.lower())
            print("Specified reference directory, reference FASTA, and protein "+\
                  "map/taxonomy table not found. " +
                  "Using database in location: " + reference_dir + ".")
            os.system("mkdir -p " + reference_dir)
            tax_tab = os.path.join(reference_dir, "tax-table.txt")
            prot_tab = os.path.join(reference_dir, "prot-map.json")
        f = open(os.path.join(output_dir, "README_DB.txt"), "a")
        e = datetime.datetime.now()
        f.write("Time finished was " + str(e) + " for database " + \
                str(args.database.lower()))

        ## Next, see whether there is a subdirectory of reference
        ## directory containing folder for our DB name
        if (not os.path.isfile(os.path.join(reference_dir, ref_fasta))) | \
           (not os.path.isfile(tax_tab)) | \
           (not os.path.isfile(prot_tab)):
            ref_fasta, tax_tab, prot_tab = downloadDatabase(args.database.lower(),
                                                            alignment_choice, output_dir,
                                                            "/".join(reference_dir.
                                                                     split("/")[0:-1]))
            if (not os.path.isfile(os.path.join(reference_dir, ref_fasta))) | \
               (not os.path.isfile(tax_tab)) | \
               (not os.path.isfile(prot_tab)):
                print("Download and formatting did not complete successfully. " +\
                      "Check log files for details.")
                sys.exit(1)
        else:
            print("Found database folder for " + reference_dir +
                  " in current directory; will not re-download.")

    if setup_choice:
        print("Creating a",alignment_choice,"reference from database files...")
        manageEukulele(piece = "setup_databases", ref_fasta = ref_fasta,
                       rerun_rules = rerun_rules, output_dir = output_dir,
                       alignment_choice = alignment_choice,
                       database_dir = reference_dir)

    if alignment_choice_select:
        ## First, we need to perform TransDecoder if needed
        manageEukulele(piece = "transdecode", mets_or_mags = mets_or_mags,
                       samples = samples, output_dir = output_dir,
                       rerun_rules = rerun_rules, sample_dir = sample_dir,
                       transdecoder_orf_size = transdecoder_orf_size,
                       nt_ext = nt_ext, pep_ext = pep_ext,
                       run_transdecoder = run_transdecoder,
                       perc_mem = perc_mem)

        ## Next to align against our database of choice ##
        alignment_res = manageEukulele(piece = "align_to_db", alignment_choice = alignment_choice,
                                       samples = samples, filter_metric = args.filter_metric,
                                       output_dir = output_dir, ref_fasta = ref_fasta,
                                       mets_or_mags = mets_or_mags, database_dir = reference_dir,
                                       sample_dir = sample_dir, rerun_rules = rerun_rules,
                                       nt_ext = nt_ext, pep_ext = pep_ext, perc_mem = perc_mem)

        ## Next to do salmon counts estimation. ##
        if use_salmon_counts:
            try:
                names_to_reads = namesToReads(reference_dir, names_to_reads, salmon_dir)
            except:
                print("The salmon directory provided could not be "+\
                      "converted to a salmon file. Check " +\
                      "above error messages, and",sys.exc_info()[0],
                      "EUKulele will continue running without counts.")
                use_salmon_counts = 0

        tax_tab_read = pd.read_csv(tax_tab, sep = "\t")
        levels_possible = ["domain","supergroup","kingdom","phylum","class","order","family","genus","species"]
        levels_file_select = list(set([curr.lower() for curr in tax_tab_read.columns \
                                          if (curr.lower() in levels_possible)]))
        levels_file = [level_curr for level_curr in levels_possible if level_curr in levels_file_select]
        manageEukulele(piece = "estimate_taxonomy", output_dir = output_dir,
                       mets_or_mags = mets_or_mags,
                       tax_tab = tax_tab, cutoff_file = args.cutoff_file,
                       consensus_cutoff = consensus_cutoff, prot_tab = prot_tab,
                       use_salmon_counts = use_salmon_counts,
                       names_to_reads = names_to_reads, alignment_res = alignment_res,
                       rerun_rules = rerun_rules, samples = samples,
                       sample_dir = sample_dir, pep_ext = pep_ext,
                       nt_ext = nt_ext, perc_mem = perc_mem,
                       level_hierarchy = levels_file)

        ## Now to visualize the taxonomy ##
        manageEukulele(piece = "visualize_taxonomy", output_dir = output_dir,
                       mets_or_mags = mets_or_mags,
                       sample_dir = sample_dir, pep_ext = pep_ext, nt_ext = nt_ext,
                       use_salmon_counts = use_salmon_counts, rerun_rules = rerun_rules,
                       level_hierarchy = levels_file)

        ## Next to assign taxonomy ##
        manageEukulele(piece = "assign_taxonomy", samples = samples, mets_or_mags = mets_or_mags,
                       sample_dir = sample_dir, pep_ext = pep_ext,
                       output_dir = output_dir,
                       level_hierarchy = levels_file)

    busco_matched = True
    if busco_choice:
        configRunBusco(output_dir = output_dir, mets_or_mags = mets_or_mags, pep_ext = pep_ext,
                       nt_ext = nt_ext, sample_dir = sample_dir, samples = samples)

        busco_matched = manageBuscoQuery(output_dir = output_dir,
                                         individual_or_summary = individual_or_summary,
                         samples = samples, mets_or_mags = mets_or_mags, pep_ext = pep_ext,
                         nt_ext = nt_ext, sample_dir = sample_dir, organisms = organisms,
                         organisms_taxonomy = organisms_taxonomy, tax_tab = tax_tab,
                         busco_threshold = args.busco_threshold, perc_mem = perc_mem)

    if core_genes & busco_matched:
        print("Investigating core genes...")
        ## Next to align against our database of choice ##
        alignment_res = manageEukulele(piece = "core_align_to_db",
                                       alignment_choice = alignment_choice,
                                       samples = samples,filter_metric = args.filter_metric,
                                       output_dir = output_dir,ref_fasta = ref_fasta,
                                       mets_or_mags = mets_or_mags, database_dir = reference_dir,
                                       sample_dir = sample_dir, rerun_rules = rerun_rules,
                                       nt_ext = nt_ext, pep_ext = pep_ext)
        if len(alignment_res) > 0:
            manageEukulele(piece = "core_estimate_taxonomy", output_dir = output_dir,
                           mets_or_mags = mets_or_mags,
                           tax_tab = tax_tab, cutoff_file = args.cutoff_file,
                           consensus_cutoff = consensus_cutoff, prot_tab = prot_tab,
                           use_salmon_counts = use_salmon_counts,
                           names_to_reads = names_to_reads, alignment_res = alignment_res,
                           rerun_rules = rerun_rules, samples = samples,
                           sample_dir = sample_dir, pep_ext = pep_ext,
                           nt_ext = nt_ext)

            ## Now to visualize the taxonomy ##
            manageEukulele(piece = "core_visualize_taxonomy", output_dir = output_dir,
                           mets_or_mags = mets_or_mags,sample_dir = sample_dir,
                           pep_ext = pep_ext, nt_ext = nt_ext,use_salmon_counts =
                           use_salmon_counts, rerun_rules = rerun_rules)

            ## Next to assign taxonomy ##
            manageEukulele(piece = "core_assign_taxonomy", samples = samples,
                           mets_or_mags = mets_or_mags,
                           sample_dir = sample_dir, pep_ext = pep_ext,
                           output_dir = output_dir)
    if not args.version:
        print("EUKulele run complete!", flush = True)

if __name__ == "__main__":
    main(args_in = " ".join(sys.argv[1:]))

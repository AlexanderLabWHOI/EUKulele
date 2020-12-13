'''
EUKulele step manager.
'''

import os
import sys
import subprocess
import multiprocessing
import math
import shutil
import pathlib
from joblib import Parallel, delayed
import pandas as pd

import EUKulele
from EUKulele.tax_placement import place_taxonomy
from EUKulele.visualize_results import visualize_all_results

from scripts.mag_stats import magStats

MEM_AVAIL_GB = 0
while MEM_AVAIL_GB == 0:
    try:
        os.system("free -m > free.csv")
        MEM_AVAIL_GB = pd.read_csv("free.csv", sep = "\s+").free[0] / 10**3
    except:
        pass

# 25 GB memory per GB file size
# add a parameter that changes between the requirement per file size (reducing to 10 GB for now)
# also add a parameter to EUKulele that decides whether you use 100% of available memory and scales
# MEM_AVAIL_GB by that amount (default to 75%)
def calc_max_jobs(num_files, size_in_bytes = 2147483648, max_mem_per_proc = 10, perc_mem = 0.75):
    size_in_gb = size_in_bytes / (1024*1024*1024)
    if size_in_gb == 0:
        size_in_gb = 0.01
    max_jobs = math.floor(MEM_AVAIL_GB * perc_mem / (max_mem_per_proc * size_in_gb * num_files)) #48)
    if max_jobs == 0:
        max_jobs = 1
    return max_jobs

max_jobs = max(1, calc_max_jobs(1))

# For DIAMOND: The program can be expected to use roughly six times this number of memory (in GB).
# So for the default value of -b2.0, the memory usage will be about 12 GB.
# So we want alignment to have -b6.0

def manageEukulele(piece, mets_or_mags = "", samples = [], database_dir = "",
                   output_dir = "", ref_fasta = "", alignment_choice = "diamond",
                   rerun_rules = False, cutoff_file = "", sample_dir = "",
                   nt_ext = "", pep_ext = "",consensus_cutoff = 0.75,
                   tax_tab = "", prot_tab = "", use_salmon_counts = False,
                   names_to_reads = "", alignment_res = "", filter_metric = "evalue",
                   run_transdecoder = False, transdecoder_orf_size = 100, perc_mem = 0.75,
                   level_hierarchy = ["domain","kingdom","phylum","class","order","family",
                                      "genus","species"]):

    """
    This function diverts management tasks to the below helper functions.
    """

    if piece == "setup_eukulele":
        setupEukulele(output_dir)
    elif piece == "setup_databases":
        createAlignmentDatabase(ref_fasta, rerun_rules, output_dir,
                                alignment_choice, database_dir)
    elif piece == "get_samples":
        return getSamples(mets_or_mags, sample_dir, nt_ext, pep_ext)
    elif piece == "transdecode":
        if mets_or_mags == "mets":
            manageTrandecode(samples, output_dir, rerun_rules, sample_dir,
                     mets_or_mags = "mets", transdecoder_orf_size = 100,
                     nt_ext = "." + nt_ext.strip('.'), pep_ext = "." + \
                             pep_ext.strip('.'),
                     run_transdecoder = run_transdecoder, perc_mem = perc_mem)
    elif piece == "align_to_db":
        return manageAlignment(alignment_choice, samples, filter_metric,
                               output_dir, ref_fasta,mets_or_mags, database_dir,
                               sample_dir, rerun_rules, nt_ext, pep_ext,
                               core = "full",perc_mem = perc_mem)
    elif piece == "estimate_taxonomy":
        manageTaxEstimation(output_dir, mets_or_mags, tax_tab, cutoff_file,
                            consensus_cutoff, prot_tab, use_salmon_counts,
                            names_to_reads, alignment_res, rerun_rules, samples,
                            sample_dir, pep_ext, nt_ext, perc_mem)
    elif piece == "visualize_taxonomy":
        manageTaxVisualization(output_dir, mets_or_mags, sample_dir, pep_ext, nt_ext,
                               use_salmon_counts, rerun_rules, level_hierarchy)
    elif piece == "assign_taxonomy":
        manageTaxAssignment(samples, mets_or_mags, output_dir, sample_dir, pep_ext, core = False)
    elif piece == "core_align_to_db":
        alignment_res = manageAlignment(alignment_choice, samples, filter_metric,
                                        output_dir, ref_fasta, mets_or_mags, database_dir,
                                        sample_dir, rerun_rules, nt_ext, pep_ext,
                                        core = "core", perc_mem = perc_mem)
        alignment_res = [curr for curr in alignment_res if curr != ""]
        return alignment_res
    elif piece == "core_estimate_taxonomy":
        manageCoreTaxEstimation(output_dir, mets_or_mags, tax_tab, cutoff_file,
                                consensus_cutoff, prot_tab, use_salmon_counts,
                                names_to_reads, alignment_res, rerun_rules, samples,
                                sample_dir, pep_ext, nt_ext, perc_mem)
    elif piece == "core_visualize_taxonomy":
        manageCoreTaxVisualization(output_dir, mets_or_mags, sample_dir, pep_ext, nt_ext,
                               use_salmon_counts, rerun_rules, level_hierarchy, core = True)
    elif piece == "core_assign_taxonomy":
        manageTaxAssignment(samples, mets_or_mags, output_dir, sample_dir, pep_ext,
                            core = True)
    else:
        print("Not a supported management function.")
        sys.exit(1)


def getSamples(mets_or_mags, sample_dir, nt_ext, pep_ext):
    """
    Get the names of the metagenomic or metatranscriptomic samples from the provided folder.
    """

    if mets_or_mags == "mets":
        samples_nt = [".".join(curr.split(".")[0:-1]) for \
                      curr in os.listdir(sample_dir) if \
                      curr.split(".")[-1] == nt_ext]
        samples_pep = [".".join(curr.split(".")[0:-1]) for \
                       curr in os.listdir(sample_dir) if \
                       curr.split(".")[-1] == pep_ext]
        samples = list(set(samples_nt + samples_pep))
        print(samples)
        if len(samples) == 0:
            print("No samples found in sample directory with",
                  "specified nucleotide or peptide extension.")
            sys.exit(1)
    else:
        samples = [".".join(curr.split(".")[0:-1]) for curr \
                   in os.listdir(sample_dir) if curr.split(".")[-1] == pep_ext]
        if len(samples) == 0:
            print("No samples found in sample directory with",
                  "specified peptide extension.")
            sys.exit(1)

    return samples


def transdecodeToPeptide(sample_name, output_dir, rerun_rules, sample_dir,
                         mets_or_mags = "mets", transdecoder_orf_size = 100,
                         nt_ext = ".fasta", pep_ext = ".faa", run_transdecoder = False):
    """
    Use TransDecoder to convert input nucleotide metatranscriptomic
    sequences to peptide sequences.
    """

    if not run_transdecoder:
        return 0
    print("Running TransDecoder for sample " + str(sample_name) + "...",
          flush = True)
    os.system("mkdir -p " + os.path.join(output_dir, mets_or_mags, "transdecoder"))
    if (os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                    sample_name + pep_ext))) & (not rerun_rules):
        print("TransDecoder file already detected for sample " +
              str(sample_name) + "; will not re-run step.", flush = True)
        return 0
    elif (os.path.isfile(os.path.join(sample_dir, sample_name + pep_ext))) & \
         (not rerun_rules):
        print("Protein files detected for sample in sample directory; " +
              "will not TransDecode.", flush = True)
        os.system("cp " + os.path.join(sample_dir, sample_name + pep_ext) + " " +
                  os.path.join(output_dir, mets_or_mags, sample_name + pep_ext))
        return 0

    TD_log = open(os.path.join(output_dir,"log","transdecoder_longorfs_" + \
                               sample_name + ".log"), "w+")
    TD_err = open(os.path.join(output_dir,"log","transdecoder_longorfs_" + \
                               sample_name + ".err"), "w+")
    if not os.path.isfile(os.path.join(sample_dir, sample_name + nt_ext)):
        print("File: " + os.path.join(sample_dir, sample_name + nt_ext) + \
              " was called by TransDecoder and "
              "does not exist. Check for typos.")
        sys.exit(1)
    rc_1 = subprocess.Popen(["TransDecoder.LongOrfs", "-t",
                            os.path.join(sample_dir, sample_name + nt_ext),
               "-m", str(transdecoder_orf_size)], stdout = TD_log,
                            stderr = TD_err).wait()
    TD_log.close()
    TD_err.close()

    TD_log = open(os.path.join(output_dir,"log","transdecoder_predict_" + \
                               sample_name + ".log"), "w+")
    TD_err = open(os.path.join(output_dir,"log","transdecoder_predict_" + \
                               sample_name + ".err"), "w+")
    rc_2 = subprocess.Popen(["TransDecoder.Predict", "-t",
                            os.path.join(sample_dir, sample_name + nt_ext),
               "--no_refine_starts"], stdout = TD_log, stderr = TD_err).wait()
    #rc_2 = p2.returncode
    TD_log.close()
    TD_err.close()

    if (rc_1 + rc_2) != 0:
        print("TransDecoder did not complete successfully for sample " +
              str(sample_name) + ". Check <output_dir>/log/ folder for details.")
        sys.exit(1)

    merged_name = sample_name + nt_ext

    os.system("mkdir -p " + os.path.join(output_dir, mets_or_mags))
    os.system("mkdir -p " + os.path.join(output_dir, mets_or_mags, "transdecoder"))

    os.replace(merged_name + ".transdecoder.pep", os.path.join(output_dir, mets_or_mags,
                                                               sample_name + pep_ext))
    os.replace(merged_name + ".transdecoder.cds", os.path.join(output_dir, mets_or_mags,
                                                               "transdecoder", sample_name +
                                                               ".fasta.transdecoder.cds"))
    os.replace(merged_name + ".transdecoder.gff3", os.path.join(output_dir, mets_or_mags,
                                                                "transdecoder", sample_name +
                                                                ".fasta.transdecoder.gff3"))
    os.replace(merged_name + ".transdecoder.bed", os.path.join(output_dir, mets_or_mags,
                                                               "transdecoder", sample_name +
                                                               ".fasta.transdecoder.bed"))
    #shutil.rmtree
    os.system("rm -rf " + merged_name + "*.transdecoder_dir*")
    return rc_1 + rc_2

def manageTrandecode(met_samples, output_dir, rerun_rules, sample_dir,
                     mets_or_mags = "mets", transdecoder_orf_size = 100,
                     nt_ext = "fasta", pep_ext = ".faa", run_transdecoder = False,
                     perc_mem = 0.75):
    """
    Now for some TransDecoding - a manager for TransDecoder steps.
    """

    if not run_transdecoder:
        return 0

    print("Running TransDecoder for MET samples...", flush = True)
    max_jobs = 1
    max_jobs_samps = [calc_max_jobs(len(met_samples),
                                    pathlib.Path(os.path.join(sample_dir,
                                                              sample + nt_ext)).stat().st_size,
                                    max_mem_per_proc = 48, perc_mem = perc_mem) \
                        for sample in met_samples  \
                        if os.path.isfile(os.path.join(sample_dir, sample + nt_ext))]
    if len(max_jobs_samps) > 0:
        max_jobs = min([calc_max_jobs(len(met_samples),
                                      pathlib.Path(os.path.join(sample_dir,
                                                                sample + nt_ext)).stat().st_size,
                                     max_mem_per_proc = 48, perc_mem = perc_mem) \
                        for sample in met_samples  \
                        if os.path.isfile(os.path.join(sample_dir, sample + nt_ext))])
    n_jobs_align = min(multiprocessing.cpu_count(), len(met_samples), max(1,max_jobs))
    transdecoder_res = Parallel(n_jobs=n_jobs_align)(delayed(transdecodeToPeptide)(sample_name,
                                                                                   output_dir,
                                                                                   rerun_rules,
                                                                                   sample_dir,
                         mets_or_mags = "mets", transdecoder_orf_size = 100,
                         nt_ext = nt_ext, pep_ext = pep_ext,
                         run_transdecoder = run_transdecoder) for sample_name in met_samples)
    all_codes = sum(transdecoder_res)
    os.system("rm -f pipeliner*")
    if all_codes > 0:
        print("TransDecoder did not complete successfully; check log folder for details.")
        sys.exit(1)
    #rcodes = [os.remove(curr) for curr in glob.glob("pipeliner*")]

def setupEukulele(output_dir):
    '''
    Perform setup tasks, such as dependency and database download.
    '''

    print("Setting things up...")
    os.system("mkdir -p " + output_dir)
    os.system("mkdir -p " + os.path.join(output_dir, "log"))

    ## Download software dependencies
    rc_1 = os.system("install_dependencies.sh references_bins/ " +
                    "1> " + os.path.join(output_dir, "log", "dependency_log.txt") + " 2> " +
                    os.path.join(output_dir, "log", "dependency_err.txt"))
    sys.path.append("references_bins/")
    os.system("echo $PATH > path_test.txt")
    if rc_1 != 0:
        print("Could not successfully install all external dependent software.\n" +
              "Check DIAMOND, BLAST, BUSCO, and TransDecoder installation.")
        return 1
    return 0

def manageAlignment(alignment_choice, samples, filter_metric, output_dir, ref_fasta,
                    mets_or_mags, database_dir, sample_dir, rerun_rules, nt_ext, pep_ext,
                    core = "full",perc_mem = 0.75):
    """
    Manage the multithreaded management of aligning to either BLAST or DIAMOND database.
    """

    print("Aligning to reference database...")
    if mets_or_mags == "mets":
        fastas = []
        for sample in samples:
            if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                           sample + "." + pep_ext)):
                fp = open(os.path.join(output_dir, mets_or_mags,
                                       sample + "." + pep_ext))
                for i, line in enumerate(fp):
                    if i == 2:
                        chars = set(list(line))
                        if len(chars) <= 5:
                            print("Peptide extension used, but this file, " +
                                  str(os.path.join(output_dir, mets_or_mags,
                                                   sample + "." + pep_ext)) +
                                  ", does not appear to be a peptide file.")
                            break
                    elif i > 2:
                        fastas.append(os.path.join(output_dir, mets_or_mags,
                                                   sample + "." + pep_ext))
                        break
                fp.close()
            elif os.path.isfile(os.path.join(sample_dir, sample + "." + pep_ext)):
                fp = open(os.path.join(sample_dir, sample + "." + pep_ext))
                for i, line in enumerate(fp):
                    if i == 2:
                        chars = set(list(line))
                        if len(chars) <= 5:
                            print("Peptide extension used, but this file, " +
                                  str(os.path.join(sample_dir, sample + "." +\
                                                   pep_ext)) +
                                  ", does not appear to be a peptide file.")
                            break
                    elif i > 2:
                        fastas.append(os.path.join(sample_dir, sample + "." +\
                                                   pep_ext))
                        break
                fastas.append(os.path.join(sample_dir, sample + "." + pep_ext))
            else:
                fastas.append(os.path.join(sample_dir, sample + "." + nt_ext))
    else:
        fastas = [os.path.join(sample_dir, sample + "." + pep_ext) for sample in samples]


    max_jobs = 1
    if len(fastas) > 0:
        max_jobs = min([calc_max_jobs(len(fastas), pathlib.Path(sample).stat().st_size,
                                     max_mem_per_proc = 10, perc_mem = perc_mem) \
                        for sample in fastas])
    n_jobs_align = min(multiprocessing.cpu_count(), len(samples), max(1,max_jobs))
    alignment_res = Parallel(n_jobs=n_jobs_align,
                             prefer="threads")(delayed(alignToDatabase)(alignment_choice,
                                                                        sample_name,
                                                                        filter_metric,
                                                                        output_dir, ref_fasta,
                                                                        mets_or_mags, database_dir,
                                                                        sample_dir, rerun_rules,
                                                                        nt_ext,pep_ext,core = core) \
                                                                    for sample_name in samples)
    #alignment_res = []
    #for sample_name in samples:
    #    alignment_res.append(alignToDatabase(alignment_choice, sample_name, filter_metric,
    #.                    output_dir, ref_fasta,
    #                     mets_or_mags, database_dir, sample_dir, rerun_rules, nt_ext, pep_ext))

    if any([((curr == None) | (curr == 1)) for curr in alignment_res]):
        print("Alignment did not complete successfully.")
        sys.exit(1)

    return alignment_res

def createAlignmentDatabase(ref_fasta, rerun_rules, output_dir,
                            alignment_choice="diamond", database_dir=""):
    """
    Creates a database from the provided reference fasta file and reference database,
    whether or not it has been autogenerated.
    """

    rc_2 = 0

    output_log = os.path.join(output_dir, "log", "alignment_out.log")
    error_log = os.path.join(output_dir, "log", "alignment_err.log")
    if alignment_choice == "diamond":
        align_db = os.path.join(database_dir, "diamond", ref_fasta.strip('.fa') + '.dmnd')
        if (not os.path.isfile(align_db)) | (rerun_rules):
            ## DIAMOND database creation ##
            os.system("mkdir -p " + os.path.join(database_dir, "diamond"))
            db = os.path.join(database_dir, "diamond", ref_fasta.strip('.fa'))
            rc_2 = os.system("diamond makedb --in " + \
                            os.path.join(database_dir, ref_fasta) + " --db " + db +\
                            " 1> " + output_log + " 2> " + error_log)
        else:
            print("Diamond database file already created; will not re-create database.",
                  flush = True)
    else:
        db = os.path.join(database_dir, "blast", ref_fasta.strip('.fa'), "database")
        os.system("mkdir -p " + db)
        db_type = "prot"
        blast_version = 5
        database = ref_fasta.strip('.')
        rc_2 = os.system("makeblastdb -in " + os.path.join(database_dir, ref_fasta) +
                        " -parse_seqids -title " + database +
                        " -dbtype " + db_type + " -out " + db + " 1> " +\
                         output_log + " 2> " +\
                         error_log)
    return rc_2

def alignToDatabase(alignment_choice, sample_name, filter_metric, output_dir, ref_fasta,
                    mets_or_mags, database_dir, sample_dir, rerun_rules,
                    nt_ext, pep_ext, core = "full"):
    """
    Align the samples against the created database.
    """

    print("Aligning sample " + sample_name + "...")
    if alignment_choice == "diamond":
        os.system("mkdir -p " + os.path.join(output_dir, mets_or_mags + "_" + core,
                                             "diamond"))
        diamond_out = os.path.join(output_dir, mets_or_mags + "_" + core,
                                   "diamond", sample_name + ".diamond.out")
        if os.path.isfile(diamond_out):
            if (pathlib.Path(diamond_out).stat().st_size != 0) & (not rerun_rules):
                print("Diamond alignment file already detected; will not re-run step.")
                return diamond_out

        align_db = os.path.join(database_dir, "diamond", ref_fasta.strip('.fa') + '.dmnd')
        alignment_method = "blastp"
        if (mets_or_mags == "mets") & (core == "full"):
            if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                           sample_name + "." + pep_ext)):
                fasta = os.path.join(output_dir, mets_or_mags,
                                     sample_name + "." + pep_ext)
            elif os.path.isfile(os.path.join(sample_dir,
                                             sample_name + "." + pep_ext)):
                fasta = os.path.join(sample_dir, sample_name + "." + pep_ext)
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)
                alignment_method = "blastx"
        elif core == "full":
            fasta = os.path.join(sample_dir, sample_name + "." + pep_ext)
        elif core == "core":
            # now concatenate the BUSCO output
            fasta = os.path.join(output_dir, sample_name + "_busco" + "." + pep_ext)
            os.system(" ".join(["concatenate_busco.sh", sample_name, fasta, output_dir]))
            if not os.path.isfile(fasta):
                print("No BUSCO matches found for sample: " + sample_name)
                return ""

        other = "--outfmt 6 -k 100 -e 1e-5"
        outfmt = 6
        k = 100
        e = 1e-5
        bitscore = 50
        pid_cutoff = 75
        diamond_log = open(os.path.join(output_dir,"log",
                                        core + "_diamond_align_" + sample_name + ".log"),
                           "w+")
        diamond_err = open(os.path.join(output_dir,"log",
                                        core + "_diamond_align_" + sample_name + ".err"),
                           "w+")
        if filter_metric == "bitscore":
            rc_1 = subprocess.Popen(["diamond", alignment_method,
                                    "--db", align_db, "-q", fasta, "-o",
                                    diamond_out, "--outfmt", str(outfmt),
                                    "-k", str(k), "--min-score",
                                   str(bitscore), '-b3.0'],
                                   stdout = diamond_log,
                                   stderr = diamond_err).wait()
            print("Diamond process exited for sample " + str(sample_name) + ".",
                  flush = True)
        elif filter_metric == "pid":
            rc_1 = subprocess.Popen(["diamond", alignment_method, "--db", align_db, "-q",
                                     fasta, "-o",diamond_out, "--outfmt", str(outfmt),
                                     "-k", str(k), "--id",str(pid_cutoff), '-b3.0'],
                                    stdout = diamond_log,
                                    stderr = diamond_err).wait()
            print("Diamond process exited for sample " + str(sample_name) + ".", flush = True)
        else:
            rc_1 = subprocess.Popen(["diamond", alignment_method,
                                     "--db", align_db, "-q", fasta, "-o",
                                     diamond_out, "--outfmt", str(outfmt),
                                     "-k", str(k), "-e", str(e), '-b3.0'],
                                    stdout = diamond_log,
                                    stderr = diamond_err).wait()
            # For debugging.
            #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #stdout, stderr = p.communicate()
            #rc_1 = p.returncode
            #print(stderr)
            #print(stdout)
            print("Diamond process exited for sample " + str(sample_name) + ".",
                  flush = True)
        if rc_1 != 0:
            print("Diamond did not complete successfully for sample",
                  str(sample_name),
                  "with rc code",str(rc_1))
            os.system("rm -f " + diamond_out)
            return 1
        return diamond_out
    else:
        blast_out = os.path.join(output_dir, mets_or_mags + "_" + core,
                                 "blast", sample_name + ".blast.txt")
        os.system("mkdir -p " + os.path.join(output_dir,
                                             mets_or_mags + "_" + core, "blast"))
        if (os.path.isfile(blast_out)) & (not rerun_rules):
            print("BLAST alignment file already detected; will not re-run step.")
            return blast_out

        align_db = os.path.join(database_dir, "blast", ref_fasta.strip('.fa'), "database")
        alignment_method = "blastp"
        if mets_or_mags == "mets":
            if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                           sample_name + "." + pep_ext)):
                fasta = os.path.join(output_dir, mets_or_mags,
                                     sample_name + "." + pep_ext)
            elif os.path.isfile(os.path.join(sample_dir, sample_name + "." + pep_ext)):
                fasta = os.path.join(sample_dir, sample_name + "." + pep_ext)
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)
                alignment_method = "blastx"
        elif core == "full":
            fasta = os.path.join(sample_dir, sample_name + "." + pep_ext)
        elif core == "core":
            # now concatenate the BUSCO output
            fasta = os.path.join(output_dir, sample_name + "_busco" + "." + pep_ext)
            os.system(" ".join(["concatenate_busco.sh", sample_name, fasta, output_dir]))
            if not os.path.isfile(fasta):
                print("No BUSCO matches found for sample: " + sample_name)
                return ""

        outfmt = 6 # tabular output format
        e = 1e-5
        os.system("export BLASTDB=" + align_db)
        blast_log = open(os.path.join(output_dir,"log",
                                      "blast_align_" + sample_name + ".log"), "w+")
        blast_err = open(os.path.join(output_dir,"log",
                                      "blast_align_" + sample_name + ".err"), "w+")
        rc_1 = subprocess.Popen([alignment_method, "-query", fasta, "-db",
                                align_db, "-out",blast_out,"-outfmt",str(outfmt),
                                "-evalue", str(e)], stdout = blast_log,
                               stderr = blast_err).wait()
        if rc_1 != 0:
            print("BLAST did not complete successfully.")
            return 1
        return blast_out

def manageTaxEstimation(output_dir, mets_or_mags, tax_tab, cutoff_file, consensus_cutoff,
                        prot_tab, use_salmon_counts, names_to_reads, alignment_res,
                        rerun_rules, samples, sample_dir, pep_ext, nt_ext, perc_mem):
    '''
    Manage the taxonomic estimation arm of the EUKulele workflow.
    '''

    print("Performing taxonomic estimation steps...", flush=True)
    os.system("mkdir -p " + os.path.join(output_dir, "taxonomy_estimation"))
    outfiles = [os.path.join(output_dir, "taxonomy_estimation", samp + \
                             "-estimated-taxonomy.out") for samp in samples]

    if mets_or_mags == "mets":
        fastas = []
        for sample in samples:
            if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                           sample + "." + pep_ext)):
                fastas.append(os.path.join(output_dir, mets_or_mags,
                                           sample + "." + pep_ext))
            elif os.path.isfile(os.path.join(sample_dir, sample + "." + pep_ext)):
                fastas.append(os.path.join(sample_dir, sample + "." + pep_ext))
            else:
                fastas.append(os.path.join(sample_dir, sample + "." + nt_ext))
    else:
        fastas = [os.path.join(sample_dir, sample + "." + pep_ext) for sample in samples]

    max_jobs = min([calc_max_jobs(len(fastas), pathlib.Path(sample).stat().st_size,
                                 max_mem_per_proc = 5, perc_mem = perc_mem) for sample in fastas])
    n_jobs_align = min(multiprocessing.cpu_count(), len(alignment_res), max(1, max_jobs))
    for t in range(len(alignment_res)):
        #curr_out = place_taxonomy(tax_tab, cutoff_file, consensus_cutoff,\
        #                                        prot_tab, use_salmon_counts, names_to_reads,\
        #                                        alignment_res[t], outfiles[t], rerun_rules)
        try:
            sys.stdout = open(os.path.join(output_dir, "log", "tax_est_" +
                                           alignment_res[t].split("/")[-1].split(".")[0] +\
                                           ".out"), "w")
            sys.stderr = open(os.path.join(output_dir, "log", "tax_est_" +
                                           alignment_res[t].split("/")[-1].split(".")[0] +\
                                           ".err"), "w")
            curr_out = place_taxonomy(tax_tab, cutoff_file, consensus_cutoff,\
                                                    prot_tab, use_salmon_counts, names_to_reads,\
                                                    alignment_res[t], outfiles[t], rerun_rules)
        except:
            print("Taxonomic estimation did not complete successfully.",
                  "Check log file for details.")
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

def manageCoreTaxEstimation(output_dir, mets_or_mags, tax_tab, cutoff_file, consensus_cutoff,
                            prot_tab, use_salmon_counts, names_to_reads, alignment_res,
                            rerun_rules, samples, sample_dir, pep_ext, nt_ext, perc_mem):
    '''
    Manage taxonomic estimation for core eukaryotic genes.
    '''

    print("Performing taxonomic estimation steps...", flush=True)
    os.system("mkdir -p " + os.path.join(output_dir, "core_taxonomy_estimation"))
    outfiles = [os.path.join(output_dir, "core_taxonomy_estimation",
                             samp + "-estimated-taxonomy.out") for samp in samples]

    if mets_or_mags == "mets":
        fastas = []
        for sample in samples:
            if os.path.isfile(os.path.join(output_dir, mets_or_mags, sample + "." + pep_ext)):
                fastas.append(os.path.join(output_dir, mets_or_mags, sample + "." + pep_ext))
            elif os.path.isfile(os.path.join(sample_dir, sample + "." + pep_ext)):
                fastas.append(os.path.join(sample_dir, sample + "." + pep_ext))
            else:
                fastas.append(os.path.join(sample_dir, sample + "." + nt_ext))
    else:
        fastas = [os.path.join(sample_dir, sample + "." + pep_ext) for sample in samples]

    max_jobs = min([calc_max_jobs(len(fastas), pathlib.Path(sample).stat().st_size,
                                 max_mem_per_proc = 10, perc_mem = perc_mem) for sample in fastas])
    n_jobs_align = min(multiprocessing.cpu_count(), len(alignment_res), max_jobs)
    for t in range(len(alignment_res)):
        try:
            sys.stdout = open(os.path.join(output_dir, "log", "core_tax_est_" +
                                           alignment_res[t].split("/")[-1].split(".")[0] + ".out"),
                              "w")
            sys.stderr = open(os.path.join(output_dir, "log", "core_tax_est_" +
                                           alignment_res[t].split("/")[-1].split(".")[0] + ".err"),
                              "w")
            curr_out = place_taxonomy(tax_tab, cutoff_file, consensus_cutoff,\
                                                    prot_tab, use_salmon_counts, names_to_reads,\
                                                    alignment_res[t], outfiles[t], rerun_rules)
        except:
            print("Taxonomic estimation for core genes did not complete successfully.",
                  "Check log file for details.")
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

def manageTaxVisualization(output_dir, mets_or_mags, sample_dir, pep_ext,\
                           nt_ext, use_salmon_counts, rerun_rules,
                           level_hierarchy):
    '''
    Management function for taxonomic visualization steps.
    '''

    print("Performing taxonomic visualization steps...", flush=True)
    out_prefix = output_dir.split("/")[-1]
    sys.stdout = open(os.path.join(output_dir, "log", "tax_vis.out"), "w")
    sys.stderr = open(os.path.join(output_dir, "log", "tax_vis.err"), "w")
    visualize_all_results(out_prefix, output_dir,
                          os.path.join(output_dir, "taxonomy_estimation"),
                          sample_dir, pep_ext, nt_ext,
                          use_salmon_counts, rerun_rules, level_hierarchy)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

def manageCoreTaxVisualization(output_dir, mets_or_mags, sample_dir,
                               pep_ext, nt_ext, use_salmon_counts,
                               rerun_rules, level_hierarchy, core = False):
    '''
    Taxonomic visualization management for core eukaryotic
    genes as identified by BUSCO.
    '''

    print("Performing taxonomic visualization steps...", flush=True)
    out_prefix = output_dir.split("/")[-1]
    sys.stdout = open(os.path.join(output_dir, "log", "core_tax_vis.out"), "w")
    sys.stderr = open(os.path.join(output_dir, "log", "core_tax_vis.err"), "w")
    visualize_all_results(out_prefix = out_prefix, out_dir = output_dir,
                          est_dir = os.path.join(output_dir, "core_taxonomy_estimation"),
                          samples_dir = sample_dir, prot_extension = pep_ext,
                          nucle_extension = nt_ext, use_counts = use_salmon_counts,
                          rerun = rerun_rules, level_hierarchy = level_hierarchy, core = core)
    #except:
    #    print("Taxonomic visualization of core genes did not complete
    #.   successfully. Check log files for details.")
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

def manageTaxAssignment(samples, mets_or_mags, output_dir, sample_dir,
                        pep_ext, core = False, perc_mem = 0.75):
    '''
    Management of taxonomic assignment process.
    '''

    if mets_or_mags == "mags":
        print("Performing taxonomic assignment steps...", flush=True)
        max_jobs = min([calc_max_jobs(len(samples),
                                      pathlib.Path(os.path.join(sample_dir,
                                                                sample + "." + \
                                                                pep_ext)).stat().st_size,
                                      max_mem_per_proc = 10,
                                      perc_mem = perc_mem)
                        for sample in samples])
        n_jobs_viz = min(multiprocessing.cpu_count(), len(samples), max(1,max_jobs))
        try:
            if core:
                assign_res = Parallel(n_jobs=n_jobs_viz,
                                      prefer="threads")(delayed(assignTaxonomy)(samp, output_dir,
                                                                                "core_taxonomy_estimation",
                                                                                mets_or_mags,
                                                                                core = True) \
                                                                           for samp in samples)
            else:
                assign_res = Parallel(n_jobs=n_jobs_viz,
                                      prefer="threads")(delayed(assignTaxonomy)(samp, output_dir,
                                                                                "taxonomy_estimation",
                                                                                mets_or_mags,
                                                                                core = False) \
                                                                           for samp in samples)
        except:
            print("Taxonomic assignment did not complete successfully.",
                  "Check log files for details.")
            sys.exit(1)

        if sum(assign_res) != 0:
            print("Taxonomic assignment did not complete successfully.",
                  "Check log files for details.")
            sys.exit(1)

def assignTaxonomy(sample_name, output_dir, est_dir, mets_or_mags, core = False):
    '''
    Assign estimated taxonomy to each MET/MAG.
    '''

    taxfile = os.path.join(output_dir, est_dir,
                           sample_name + "-estimated-taxonomy.out")
    levels_directory = os.path.join(output_dir,
                                    "levels_mags")
    max_dir = os.path.join(output_dir, "max_level_mags")
    if core:
        levels_directory = os.path.join(output_dir,
                                        "core_levels_mags")
        max_dir = os.path.join(output_dir,
                               "core_max_level_mags")

    error_log = os.path.join(output_dir, "log",
                             "_".join(est_dir.split("_")[0:1]) +\
                             "_assign_" + sample_name + ".err")
    out_log = os.path.join(output_dir, "log",
                           "_".join(est_dir.split("_")[0:1]) +\
                           "_assign_" + sample_name + ".out")

    sys.stdout = open(out_log, "w")
    sys.stderr = open(error_log, "w")
    try:
        rc = magStats(["--estimated-taxonomy-file",taxfile,
                       "--out-prefix",sample_name,"--outdir",
                       levels_directory,"--max-out-dir",max_dir])
    except:
        print("Taxonomic assignment did not complete successfully for sample",
              str(sample_name),". Check log for details.")
        sys.exit(1)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return rc

'''
Facilitates the call to the BUSCO software from within EUKulele.
'''

import os
import sys
import math
import subprocess
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd

from scripts.query_busco import queryBusco

MEM_AVAIL_GB = 0
while MEM_AVAIL_GB == 0:
    try:
        os.system("free -m > free.csv")
        MEM_AVAIL_GB = pd.read_csv("free.csv", sep = "\s+").free[0] / 10**3
    except:
        pass

def calc_max_jobs(num_files, size_in_bytes = 2147483648, max_mem_per_proc = 40, perc_mem = 0.75):
    '''
    Calculate the maximum number of jobs to be run simultaneously on the system.
    '''

    size_in_gb = size_in_bytes / (1024*1024*1024)
    if size_in_gb == 0:
        size_in_gb = 0.01
    max_jobs = math.floor(MEM_AVAIL_GB * perc_mem / \
                          (max_mem_per_proc * size_in_gb * num_files))
    if max_jobs == 0:
        max_jobs = 1
    return max_jobs

max_jobs = calc_max_jobs(10)

def readBuscoFile(individual_or_summary, busco_file, organisms, organisms_taxonomy):
    if individual_or_summary == "individual":
        if (busco_file != "") & (os.path.isfile(busco_file)):
            busco_file_read = pd.read_csv(busco_file, sep = "\t")
            organisms = list(busco_file_read.iloc[:,0])
            organisms_taxonomy = list(busco_file_read.iloc[:,1])
            print("Organisms and their taxonomy levels for BUSCO analysis were read from file.")

            if len(organisms) != len(organisms_taxonomy):
                print("Organisms and taxonomic specifications for BUSCO analysis "+\
                      "do not contain the same number of entries. " + \
                      "Please revise such that each organism flagged for "+\
                      "BUSCO analysis also includes its original " + \
                      "taxonomic level.")
                sys.exit(1)
        else:
            print("No BUSCO file specified/found; using argument-specified " +\
                  "organisms and taxonomy for BUSCO analysis.")

    return organisms, organisms_taxonomy


def configRunBusco(output_dir, mets_or_mags, pep_ext, nt_ext, sample_dir, samples):
    print("Performing BUSCO steps...", flush=True)
    print("Configuring BUSCO...", flush=True)

    ## Run BUSCO on the full dataset ##
    busco_db = "eukaryota_odb10"
    busco_config_res = configure_busco(busco_db,output_dir)
    n_jobs_busco = min(multiprocessing.cpu_count(), len(samples),
                       max(1, calc_max_jobs(len(samples))))
    print("Running busco with",n_jobs_busco,"simultaneous jobs...", flush=True)
    busco_res = Parallel(n_jobs=n_jobs_busco, prefer="threads")\
                        (delayed(run_busco)(sample_name,os.path.join(output_dir,"busco"),
                                            output_dir,busco_db, mets_or_mags,pep_ext,
                                            nt_ext,sample_dir) \
                         for sample_name in samples)
    print(os.listdir(os.path.join(output_dir, "busco", samples[0])),
          "is what is in BUSCO directory")
    all_codes = sum(busco_res) + busco_config_res
    if sum(busco_res) > 0:
        print("BUSCO initial run did not complete successfully.\n" +
              "Please check the BUSCO run log files in the log/ folder.", flush = True)
        sys.exit(1)
    if busco_config_res > 0:
        print("BUSCO initial configuration did not complete successfully.\n" +
              "Please check the BUSCO configuration log files in the log/ folder.", flush = True)
        sys.exit(1)

def configure_busco(busco_db,output_dir):
    busco_config_log = open(os.path.join(output_dir,"log",
                                         "busco_config.out"), "w+")
    busco_config_err = open(os.path.join(output_dir,"log",
                                         "busco_config.err"), "w+")
    rc1 = 0

    if not os.path.isdir(os.path.join("busco_downloads","lineages",
                                      "eukaryota_odb10")):
        p1 = subprocess.Popen(["configure_busco.sh", busco_db],
                              stdout = busco_config_log, stderr = busco_config_err)
        p1.wait()
        rc1 = p1.returncode
    else:
        print("BUSCO lineage database already found; not re-downloaded.")

    busco_config_log.close()
    busco_config_err.close()
    return rc1

def run_busco(sample_name, output_dir_busco, output_dir, busco_db,
              mets_or_mags, pep_ext, nt_ext, sample_dir):
    '''
    Runs the BUSCO software.
    '''

    CPUS = multiprocessing.cpu_count()

    if mets_or_mags == "mets":
        if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                       sample_name + "." + pep_ext)):
            fastaname = os.path.join(output_dir, mets_or_mags,
                                     sample_name + "." + pep_ext)
            busco_mode = "proteins"
        elif os.path.isfile(os.path.join(sample_dir, sample_name + "." + pep_ext)):
            fastaname = os.path.join(sample_dir, sample_name + "." + pep_ext)
            busco_mode = "proteins"
        else:
            fastaname = os.path.join(sample_dir, sample_name + "." + nt_ext)
            busco_mode = "transcriptome"
    else:
        fastaname = os.path.join(sample_dir, sample_name + "." + pep_ext)
        busco_mode = "proteins"

    busco_run_log = open(os.path.join(output_dir,"log","busco_run.out"), "w+")
    busco_run_err = open(os.path.join(output_dir,"log","busco_run.err"), "w+")
    p1 = subprocess.Popen(["run_busco.sh", str(sample_name), str(output_dir_busco),
                              os.path.join(output_dir_busco, "config_" +
                                           sample_name + ".ini"),
                              fastaname, str(CPUS), busco_db, busco_mode],
                          stdout = busco_run_log, stderr = busco_run_err)

    ## TRAVIS DEBUGGING!! ##

    p1.wait()
    rc1 = p1.returncode

    busco_run_log.close()
    busco_run_err.close()


    #a_file = open(os.path.join(output_dir,"log","busco_run.err"))

    #lines = a_file.readlines()
    #print("BUSCO error log:")
    #for line in lines:
    #    print(line)

    #a_file = open(os.path.join(output_dir,"log","busco_run.out"))

    #lines = a_file.readlines()
    #print("BUSCO output log:")
    #for line in lines:
    #    print(line)

    return rc1

def manageBuscoQuery(output_dir, individual_or_summary, samples,
                     mets_or_mags, pep_ext, nt_ext,
                     sample_dir, organisms, organisms_taxonomy,
                     tax_tab, busco_threshold, perc_mem):
    """
    Assess BUSCO completeness on the most prevalent members of the
    metatranscriptome at each taxonomic level.
    """
    max_jobs = calc_max_jobs(len(samples), perc_mem = perc_mem)
    samples_complete = []
    if individual_or_summary == "individual":
        if len(organisms) != len(organisms_taxonomy):
            print("A different number of organisms was specified than "+\
                  "the taxonomic levels given in " +\
                  "individual mode. Please check inputs.")
            sys.exit(1)
        if (len(organisms) == 0) | (len(organisms_taxonomy) == 0):
            print("The number of organisms specified was " + str(len(organisms)) +\
                  " and the number of taxonomic levels specified was " +\
                  str(len(organisms_taxonomy)) +
                  " in individual mode. Neither can be zero. "+\
                  "Please check inputs.")
            sys.exit(1)

        for sample_name in samples:
            # the BUSCO table that we're interested in using that contains the
            # BUSCO matches and their level of completeness
            if not os.path.isfile(os.path.join(output_dir, "busco",
                                               sample_name, "run_eukaryota_odb10",
                                               "full_table.tsv")):
                print("BUSCO run either did not complete successfully, "+\
                      "or returned no matches for sample",
                      sample_name,". Check busco_run log for details.")
                continue
            samples_complete.append(sample_name)
 
            busco_table = os.path.join(output_dir, "busco", sample_name,
                                       "run_eukaryota_odb10", "full_table.tsv")
            missing_buscos = pd.read_csv(os.path.join(output_dir, "busco",
                                                      sample_name,
                                                      "run_eukaryota_odb10",
                                                      "missing_busco_list.tsv"),
                                         sep = "\t", comment = "#", header = None)
            if len(missing_buscos.index) < 255:
                print("At least one BUSCO present in sample",sample_name,"but",
                      len(missing_buscos.index),
                      "missing.",flush=True)
                samples_complete.append(sample_name)
            else:
                print("No matches returned for sample",sample_name,
                      ". Assessment files will be empty.",flush=True)
            # the prefix to specify where the taxonomy estimation
            # output files are located
            taxfile_stub = os.path.join(output_dir, "taxonomy_counts",
                                        output_dir.split("/")[-1])

            if mets_or_mags == "mets":
                if os.path.isfile(os.path.join(output_dir, mets_or_mags,
                                               sample_name + "." + pep_ext)):
                    fasta = os.path.join(output_dir, mets_or_mags,
                                         sample_name + "." + pep_ext)
                else:
                    fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + pep_ext)

            query_busco_log = open(os.path.join(output_dir,"log","busco_query_" +\
                                                sample_name + ".log"), "w+")
            query_busco_err = open(os.path.join(output_dir,"log","busco_query_" +\
                                                sample_name + ".err"), "w+")
            sys.stdout = query_busco_log
            sys.stderr = query_busco_err
            query_args = ["--organism_group",str(" ".join(organisms)),"--taxonomic_level",
                          str(" ".join(organisms_taxonomy)),"--output_dir",output_dir,
                          "--fasta_file",fasta,"--sample_name",sample_name,
                          "--taxonomy_file_prefix",taxfile_stub,
                          "--tax_table",tax_tab,"--busco_out",busco_table,
                          "-i","individual",
                          "--busco_threshold",str(busco_threshold)]
            try:
                rc = queryBusco(query_args)
            except:
                print("BUSCO query did not run successfully for sample " +\
                      sample_name + "; check log file for details.")
                sys.exit(1)

            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            #if (rc != 0) | (os.stat(os.path.join("log","busco_query_" + sample_name +
            # ".err")).st_size != 0):
            #    print("BUSCO query did not run successfully for sample
            # " + sample_name + "; check log file for details.")
            #    sys.exit(1)
    else:
        for sample_name in samples:
            # the BUSCO table that we're interested in using that contains the
            # BUSCO matches and their level of completeness
            if not os.path.isfile(os.path.join(output_dir, "busco", sample_name,
                                               "run_eukaryota_odb10", "full_table.tsv")):
                print("BUSCO run either did not complete successfully, ",
                      "or returned no matches for sample",
                      sample_name,". Check busco_run log for details.",
                      flush=True)
                continue
            busco_table = os.path.join(output_dir, "busco", sample_name,
                                       "run_eukaryota_odb10", "full_table.tsv")
            missing_buscos = pd.read_csv(os.path.join(output_dir, "busco", sample_name,
                                                      "run_eukaryota_odb10",
                                                      "missing_busco_list.tsv"),
                                         sep = "\t", comment = "#", header = None)
            if len(missing_buscos.index) < 255:
                print("At least one BUSCO present in sample",sample_name,
                      "but",len(missing_buscos.index),
                      "missing.",flush=True)
                samples_complete.append(sample_name)
            else:
                print("No matches returned for sample",sample_name,
                      ". Assessment files will be empty.",flush=True)
            # the prefix to specify where the taxonomy estimation output files are located
            taxfile_stub = os.path.join(output_dir, "taxonomy_counts",
                                        output_dir.split("/")[-1])

            if mets_or_mags == "mets":
                fasta = os.path.join(output_dir, mets_or_mags,
                                     sample_name + "." + pep_ext)
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)

            query_busco_log = open(os.path.join(output_dir,"log",
                                                "busco_query_" + sample_name + ".log"),
                                   "w+")
            query_busco_err = open(os.path.join(output_dir,"log",
                                                "busco_query_" + sample_name + ".err"),
                                   "w+")
            sys.stdout = query_busco_log
            sys.stderr = query_busco_err
            query_args = ["--output_dir",output_dir,"--fasta_file",fasta,"--sample_name",
                          sample_name,"--taxonomy_file_prefix",taxfile_stub,"--tax_table",
                          tax_tab,"--busco_out",busco_table,"-i","summary"]

            try:
                rc = queryBusco(query_args)
            except OSError as e:
                print("Not all files needed to run BUSCO query ",
                      "(output of BUSCO run) found;",\
                      "check log file for details. Here is the error:",e)
                rc = 1
            except:
                print("Unexpected error:", sys.exc_info()[0])
                rc = 1

            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

            if (rc != 0) | (os.stat(os.path.join(output_dir,"log",
                                                 "busco_query_" +\
                                                 sample_name + ".err")).st_size != 0):
                print("BUSCO query did not run successfully for sample " +\
                      sample_name + "; check log file for details.")
                sys.exit(1)
            else:
                print("BUSCO query complete.")
   
    if len(samples_complete) == 0:
        print("No BUSCO matches found for any sample. ",
              "Check BUSCO run log for details. Exiting...")
        return False
    return True

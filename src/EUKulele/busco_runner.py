import os
import sys
import subprocess
import multiprocessing
from joblib import Parallel, delayed

from scripts.query_busco import queryBusco

def readBuscoFile(individual_or_summary, busco_file, organisms, organisms_taxonomy):
    if individual_or_summary == "individual":
        if (busco_file != "") & (os.path.isfile(busco_file)):
            busco_file_read = read.csv(busco_file, sep = "\t")
            organisms = list(busco_file_read.iloc[:,0])
            organisms_taxonomy = list(busco_file_read.iloc[:,1])
            print("Organisms and their taxonomy levels for BUSCO analysis were read from file.")
            
            if (len(organisms) != len(organisms_taxonomy)):
                print("Organisms and taxonomic specifications for BUSCO analysis do not contain the same number of entries. " + \
                      "Please revise such that each organism flagged for BUSCO analysis also includes its original " + \
                      "taxonomic level.")
                sys.exit(1)
        else:
            print("No BUSCO file specified/found; using argument-specified organisms and taxonomy for BUSCO analysis.")
            
    return organisms, organisms_taxonomy


def configRunBusco(output_dir, mets_or_mags, pep_ext, nt_ext, sample_dir, samples):
    print("Performing BUSCO steps...", flush=True)
    print("Running busco...")
    
    ## Run BUSCO on the full dataset ##
    busco_db = "eukaryota_odb10"
    busco_config_res = configure_busco(busco_db)
    n_jobs_busco = min(multiprocessing.cpu_count(), len(samples), 4)
    busco_res = Parallel(n_jobs=n_jobs_busco, prefer="threads")(delayed(run_busco)(sample_name, 
                                                                                                  os.path.join(output_dir, 
                                                                                                               "busco"), 
                                                                                                  output_dir,
                                                                                                  busco_db, mets_or_mags, 
                                                                                                  pep_ext, nt_ext,
                                                                                                  sample_dir) \
                                                                               for sample_name in samples)
    all_codes = sum(busco_res) + busco_config_res
    if sum(busco_res) > 0:
        print("BUSCO initial run did not complete successfully.\n" + 
              "Please check the BUSCO run log files in the log/ folder.")
        sys.exit(1)
    if busco_config_res > 0:
        print("BUSCO initial configuration did not complete successfully.\n" + 
              "Please check the BUSCO configuration log files in the log/ folder.")
        sys.exit(1)
                
def configure_busco(busco_db):
    busco_config_log = open(os.path.join("log","busco_config.out"), "w+")
    busco_config_err = open(os.path.join("log","busco_config.err"), "w+")
    rc1 = 0
    
    if not os.path.isdir(os.path.join("busco_downloads","lineages","eukaryota_odb10")):
        p1 = subprocess.Popen(["configure_busco.sh", busco_db], stdout = busco_config_log, stderr = busco_config_err)
        p1.wait()
        rc1 = p1.returncode
    else:
        print("BUSCO lineage database already found; not re-downloaded.")
        
    busco_config_log.close()
    busco_config_err.close()
    return rc1
        
def run_busco(sample_name, output_dir_busco, output_dir, busco_db, mets_or_mags, pep_ext, nt_ext, sample_dir):
    CPUS = multiprocessing.cpu_count()
    if mets_or_mags == "mets":
        fastaname = os.path.join(output_dir, mets_or_mags, sample_name + "." + pep_ext) 
    else:
        fastaname = os.path.join(sample_dir, sample_name + "." + pep_ext)
        
    busco_run_log = open(os.path.join("log","busco_run.out"), "w+")
    busco_run_err = open(os.path.join("log","busco_run.err"), "w+")
    p1 = subprocess.Popen(["run_busco.sh", str(sample_name), str(output_dir_busco), 
                              os.path.join(output_dir_busco, "config_" + sample_name + ".ini"), 
                              fastaname, str(CPUS), busco_db], stdout = busco_run_log, stderr = busco_run_err)
    p1.wait()
    rc1 = p1.returncode
    
    busco_run_log.close()
    busco_run_err.close()
    return rc1 

def manageBuscoQuery(output_dir, individual_or_summary, samples, mets_or_mags, pep_ext, nt_ext,
                     sample_dir, organisms, organisms_taxonomy, tax_tab):
    """
    Assess BUSCO completeness on the most prevalent members of the metatranscriptome at each taxonomic level.
    """
    
    if individual_or_summary == "individual":
        for sample_name in samples:
            # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
            busco_table = os.path.join(output_dir, "busco", sample_name, "full_table.tsv") 
            # the prefix to specify where the taxonomy estimation output files are located
            taxfile_stub = os.path.join(output_dir, "taxonomy_counts", output_dir.split("/")[-1]) 

            if mets_or_mags == "mets":
                fasta = os.path.join(output_dir, mets_or_mags, sample_name + "." + pep_ext) 
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)

            query_busco_log = open(os.path.join("log","busco_query_" + sample_name + ".log"), "w+")
            query_busco_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
            sys.stdout = query_busco_log
            sys.stderr = query_busco_err
            query_args = ["--organism_group",str(" ".join(organisms)),"--taxonomic_level",
                          str(" ".join(organisms_taxonomy)),"--output_dir",output_dir,"--fasta_file",
                          fasta,"--sample_name",sample_name,"--taxonomy_file_prefix",taxfile_stub,
                          "--tax_table",tax_tab,"--busco_out",busco_table,"-i","individual"]
            rc = queryBusco(query_args)

            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__ 
            if (rc != 0) | (os.stat(os.path.join("log","busco_query_" + sample_name + ".err")).st_size != 0):
                print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")
                sys.exit(1)
    else:
        for sample_name in samples:
            # the BUSCO table that we're interested in using that contains the BUSCO matches and their level of completeness
            busco_table = os.path.join(output_dir, "busco", sample_name, "full_table.tsv") 
            # the prefix to specify where the taxonomy estimation output files are located
            taxfile_stub = os.path.join(output_dir, "taxonomy_counts", output_dir.split("/")[-1])

            if mets_or_mags == "mets":
                fasta = os.path.join(output_dir, mets_or_mags, sample_name + "." + pep_ext) 
            else:
                fasta = os.path.join(sample_dir, sample_name + "." + nt_ext)

            query_busco_log = open(os.path.join("log","busco_query_" + sample_name + ".log"), "w+")
            query_busco_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
            sys.stdout = query_busco_log
            sys.stderr = query_busco_err
            query_args = ["--output_dir",output_dir,"--fasta_file",fasta,"--sample_name",
                          sample_name,"--taxonomy_file_prefix",taxfile_stub,"--tax_table",
                          tax_tab,"--busco_out",busco_table,"-i","summary"]

            try:
                rc = queryBusco(query_args)
            except OSError as e:
                print("Not all files needed to run BUSCO query (output of BUSCO run) not found; check log file for details.")
            except:
                print("Unexpected error:", sys.exc_info()[0])
                
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__ 

            if (rc != 0) | (os.stat(os.path.join("log","busco_query_" + sample_name + ".err")).st_size != 0):
                print("BUSCO query did not run successfully for sample " + sample_name + "; check log file for details.")
                sys.exit(1)
            else:
                print("BUSCO query complete.")
import os
import sys
import subprocess

def manageEukulele(piece, mets_or_mags = "", samples = [], database_dir = "", 
                   output_dir = "", ref_fasta = "", alignment_choice = "diamond", 
                   database_dir = "", rerun_rules = False):
    """
    This function diverts management tasks to the below helper functions.
    """
    
    if piece == "setup_eukulele":
        setupEukulele()
    elif piece == "setup_databases":
        createAlignmentDatabase(ref_fasta, rerun_rules, alignment_choice, database_dir)
    elif piece == "get_samples":
        return getSamples(mets_or_mags)
    elif piece == "transdecode":
        if mets_or_mags == "mets":
            manageTrandecode(samples)
    elif piece == "align_to_db":
        manageAlignment(alignment_choice, sample_names, filter_metric, output_dir, ref_fasta, 
                        mets_or_mags, database_dir, sample_dir, rerun_rules, nt_ext, pep_ext)
        
             
def getSamples():
    """
    Get the names of the metagenomic or metatranscriptomic samples from the provided folder.
    """
    
    if (mets_or_mags == "mets"):
        samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == NT_EXT]
        if len(samples) == 0:
            print("No samples found in sample directory with specified nucleotide extension.")
            sys.exit(1)
    else:
        samples = [".".join(curr.split(".")[0:-1]) for curr in os.listdir(SAMPLE_DIR) if curr.split(".")[-1] == PEP_EXT]
        if len(samples) == 0:
            print("No samples found in sample directory with specified peptide extension.")
            sys.exit(1)
            
    return samples
            

def transdecodeToPeptide(sample_name):
    """
    Use TransDecoder to convert input nucleotide metatranscriptomic sequences to peptide sequences.
    """
    
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder"))
    if (os.path.isfile(os.path.join(OUTPUTDIR, mets_or_mags, 
                                    sample_name + "." + PEP_EXT))) & (not RERUN_RULES):
        print("TransDecoder file already detected for sample " + 
              str(sample_name) + "; will not re-run step.")
        return 0
    
    TD_log = open(os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT), "w+")
    TD_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
    p1 = Popen(["TransDecoder.LongOrfs", "-t", os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT),
               "-m", str(TRANSDECODERORFSIZE)], stdout = TD_log, stderr = TD_err)
    p1.wait()
    rc1 = p1.returncode
    
    TD_log = open(os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT), "w+")
    TD_err = open(os.path.join("log","busco_query_" + sample_name + ".err"), "w+")
    p2 = Popen(["TransDecoder.Predict -t", "-t", os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT),
               "--no_refine_starts"], stdout = TD_log, stderr = TD_err)
    p2.wait()
    rc2 = p2.returncode
    
    if (rc1 + rc2) != 0:
        print("TransDecoder did not complete successfully for sample " + 
              str(sample_name) + ". Check log/ folder for details.")
        sys.exit(1)
        
    merged_name = sample_name + "." + NT_EXT
    
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags))
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder"))
    
    os.replace(merged_name + ".transdecoder.pep", os.path.join(OUTPUTDIR, mets_or_mags, 
                                                               sample_name + "." + PEP_EXT))
    os.replace(merged_name + ".transdecoder.cds", os.path.join(OUTPUTDIR, mets_or_mags, 
                                                               "transdecoder", sample_name + 
                                                               ".fasta.transdecoder.cds"))
    os.replace(merged_name + ".transdecoder.gff3", os.path.join(OUTPUTDIR, mets_or_mags, 
                                                                "transdecoder", sample_name + 
                                                                ".fasta.transdecoder.gff3"))
    os.replace(merged_name + ".transdecoder.bed", os.path.join(OUTPUTDIR, mets_or_mags, 
                                                               "transdecoder", sample_name + 
                                                               ".fasta.transdecoder.bed"))
    shutil.rmtree(merged_name + ".transdecoder_dir*")
    return rc1 + rc2
    
def manageTrandecode(met_samples):
    """
    Now for some TransDecoding - a manager for TransDecoder steps.
    """
    
    n_jobs_align = min(multiprocessing.cpu_count(), len(met_samples))
    transdecoder_res = Parallel(n_jobs=n_jobs_align)(delayed(transdecodeToPeptide)(sample_name) for sample_name in met_samples)
    all_codes = sum(transdecoder_res)
    if all_codes > 0:
        print("TransDecoder did not complete successfully; check log folder for details.")
        sys.exit(1)
    rcodes = [os.remove(curr) for curr in glob.glob("pipeliner*")]
              
def setupEukulele():
    print("Setting things up...")
    os.system("mkdir -p " + OUTPUTDIR)
    os.system("mkdir -p log")

    ## Download software dependencies
    rc1 = os.system("source " + "install_dependencies.sh references_bins/ " + 
                    "1> log/dependency_log.txt 2> log/dependency_err.txt")
    sys.path.append("references_bins/")
    os.system("echo $PATH > path_test.txt")
    if rc1 != 0:
        print("Could not successfully install all external dependent software.\n" + 
              "Check DIAMOND, BLAST, BUSCO, and TransDecoder installation.")
        return 1
    return 0

def manageAlignment(alignment_choice, sample_name, filter_metric, output_dir, ref_fasta,
                    mets_or_mags, database_dir, sample_dir, rerun_rules, nt_ext, pep_ext):
    """
    Manage the multithreaded management of aligning to either BLAST or DIAMOND database.
    """
    
    n_jobs_align = min(multiprocessing.cpu_count(), len(samples))
    alignment_res = Parallel(n_jobs=n_jobs_align, prefer="threads")(delayed(align_to_database)(alignment_choice,
                                                                                               sample_name, filter_metric, 
                                                                                               output_dir, ref_fasta, 
                                                                                               mets_or_mags, database_dir, 
                                                                                               sample_dir, rerun_rules, nt_ext, 
                                                                                               pep_ext) \
                                                                    for sample_name in samples)
    
    if any([((curr == None) | (curr == 1)) for curr in alignment_res]):
        print("Alignment did not complete successfully.")
        sys.exit(1)

def createAlignmentDatabase(REF_FASTA, RERUN_RULES, alignment_choice="diamond", DATABASE_DIR=""):
    """
    Creates a database from the provided reference fasta file and reference database,
    whether or not it has been autogenerated.
    """
    
    rc2 = 0
              
    output_log = "alignment_out.log"
    error_log = "alignment_err.log"
    if alignment_choice == "diamond":
        align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
        if (not os.path.isfile(align_db)) | (RERUN_RULES):
            ## DIAMOND database creation ##
            db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa'))
            rc2 = os.system("diamond makedb --in " + REF_FASTA + " --db " + db + 
                            " 1> " + output_log + " 2> " + error_log)
        else:
            print("Diamond database file already created; will not re-create database.", flush = True)
    else:
        db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
        db_type = "prot"
        blast_version = 5
        rc2 = os.system("makeblastdb -in " + REF_FASTA + " -parse_seqids -blastdb_version " + 
                        str(blast_version) + " -title " + database + " -dbtype " + db_type + " -out " + db)
    return rc2
 
def alignToDatabase(alignment_choice, sample_name, filter_metric, OUTPUTDIR, REF_FASTA, 
                      mets_or_mags, DATABASE_DIR, SAMPLE_DIR, RERUN_RULES, NT_EXT, PEP_EXT):
    """
    Align the 
    """
    
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
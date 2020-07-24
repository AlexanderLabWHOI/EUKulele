import os
    
def setupEukulele():
    print("Setting things up...")
    os.system("mkdir -p " + OUTPUTDIR)
    os.system("mkdir -p log")

    ## Download software dependencies
    rc1 = os.system("source " + "install_dependencies.sh references_bins/ 1> log/dependency_log.txt 2> log/dependency_err.txt")
    sys.path.append("references_bins/")
    os.system("echo $PATH > path_test.txt")
    if rc1 != 0:
        print("Could not successfully install all external dependent software. Check DIAMOND, BLAST, BUSCO, and TransDecoder installation.")
        return 1
    return 0


def setupDatabases(REF_FASTA, RERUN_RULES, alignment_choice="diamond", DATABASE_DIR=""):
    rc2 = 0
    if DATABASE_DIR != "":
        output_log = "alignment_out.log"
        error_log = "alignment_err.log"
        if alignment_choice == "diamond":
            align_db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa') + '.dmnd')
            if (not os.path.isfile(align_db)) | (RERUN_RULES):
                ## DIAMOND database creation ##
                db = os.path.join(DATABASE_DIR, "diamond", REF_FASTA.strip('.fa'))
                rc2 = os.system("diamond makedb --in " + concatenated_file + " --db " + db + " 1> " + output_log + " 2> " + error_log)
            else:
                print("Diamond database file already created; will not re-create database.", flush = True)
        else:
            db = os.path.join(DATABASE_DIR, "blast", REF_FASTA.strip('.fa'), "database")
            db_type = "prot"
            blast_version = 5
            rc2 = os.system("makeblastdb -in " + concatenated_file + " -parse_seqids -blastdb_version " + str(blast_version) + " -title " + database + " -dbtype " + db_type + " -out " + db)
    return rc2
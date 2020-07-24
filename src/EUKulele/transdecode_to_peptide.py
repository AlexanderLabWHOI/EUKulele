import os

def transdecodeToPeptide(sample_name):
    ## THIS ALL NEEDS TO BE A CALLED BASH SCRIPT - look into subprocess.Popen with bash arguments
    os.system("mkdir -p " + os.path.join(OUTPUTDIR, mets_or_mags, "transdecoder"))
    if (os.path.isfile(os.path.join(OUTPUTDIR, mets_or_mags, sample_name + "." + PEP_EXT))) & (not RERUN_RULES):
        print("TransDecoder file already detected for sample " + str(sample_name) + "; will not re-run step.")
        return 0
    rc1 = os.system("TransDecoder.LongOrfs -t " + os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT) + " -m " + str(TRANSDECODERORFSIZE) + " 2> " + os.path.join("log", "transdecoder_error_" + sample_name + ".err") + " 1> " + os.path.join("log", "transdecoder_out_" + sample_name + ".out"))
    rc2 = os.system("TransDecoder.Predict -t " + os.path.join(SAMPLE_DIR, sample_name + "." + NT_EXT) + " --no_refine_starts 2>> " + os.path.join("log", "transdecoder_error_" + sample_name + ".err") + " 1>> " + os.path.join("log", "transdecoder_out_" + sample_name + ".out"))
    
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
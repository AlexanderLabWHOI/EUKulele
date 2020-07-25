import os
import sys

def Manager(piece, mets_or_mags = "", ):
    if piece == "get_samples":
        return getSamples(mets_or_mags)
    elif piece == "transdecode":
        if mets_or_mags == "mets":
            manageTrandecode(samples)
    elif piece == 
     
        
def getSamples():
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
    
def manageTrandecode(met_samples):
    ## Now for some TransDecoding ##
    n_jobs_align = min(multiprocessing.cpu_count(), len(met_samples))
    transdecoder_res = Parallel(n_jobs=n_jobs_align)(delayed(transdecodeToPeptide)(sample_name) for sample_name in met_samples)
    all_codes = sum(transdecoder_res)
    if all_codes > 0:
        print("TransDecoder did not complete successfully; check log folder for details.")
        sys.exit(1)
    rcodes = [os.remove(curr) for curr in glob.glob("pipeliner*")]
 
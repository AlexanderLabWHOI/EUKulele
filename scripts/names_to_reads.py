''' Function for dealing with Salmon imput to EUKulele.'''

import os
import sys
import glob
import pandas as pd

def namesToReads(reference_dir, names_to_reads, salmon_dir):
    ''' Main function to create reads table from Salmon names.'''

    if os.path.isfile(os.path.join(reference_dir,names_to_reads)):
        print("Salmon reads file previously created; new file will" +\
              "not be created from Salmon directory.")
        sys.exit(0)

    folder_names = glob.glob(os.path.join(salmon_dir,'*quant*'))
    files_salmon = [os.path.join(curr,"quant.sf") for curr in folder_names]

    transcript_dict = dict()
    transcript_sample_dict = dict()

    for curr_ind in range(len(folder_names)):
        curr_salmon = pd.read_csv(files_salmon[curr_ind], sep = "\t")
        for curr_ind_2 in range(len(curr_salmon.index)):
            name_curr = curr_salmon["Name"][curr_ind_2]
            read_curr = float(curr_salmon["NumReads"][curr_ind_2])
            sample_curr = folder_names[curr_ind].split("_")[-1]
            if name_curr in transcript_dict:
                transcript_dict[name_curr] = transcript_dict[name_curr] + read_curr
                transcript_sample_dict[name_curr].append(sample_curr)
            else:
                transcript_dict[name_curr] = read_curr
                transcript_sample_dict[name_curr] = [sample_curr]

    names_to_reads = pd.DataFrame({"TranscriptNames": list(transcript_dict.keys()),
                                   "NumReads": list(transcript_dict.values()),
                                   "SampleName": list(transcript_sample_dict.values())})

    if ".csv" in names_to_reads:
        names_to_reads.to_csv(path_or_buf =
                              os.path.join(reference_dir,names_to_reads), sep = "\t")
    else:
        names_to_reads.to_csv(path_or_buf =
                              os.path.join(reference_dir,"namestoreads.csv"), sep = "\t")
        names_to_reads = os.path.join(reference_dir,"namestoreads.csv")

    return names_to_reads

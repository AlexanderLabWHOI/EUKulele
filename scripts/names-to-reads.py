import pandas as pd
import glob
import os
import yaml
import sys

with open("config.yaml", 'r') as configfile:
    config = yaml.safe_load(configfile)

if os.path.isfile(os.path.join(config["reference"],config["names_to_reads"])):
    print("Salmon reads file previously created; new file will not be created from Salmon directory.")
    sys.exit(0)

salmon_dir = config["salmon_dir"]
folder_names = glob.glob(os.path.join(salmon_dir,'*quant*'))
files_salmon = [os.path.join(curr,"quant.sf") for curr in folder_names]

transcript_names = []
transcript_counts = []
sample_names = []

for curr_ind in range(len(folder_names)):
    curr_salmon = pd.read_csv(files_salmon[curr_ind], sep = "\t")
    transcript_names.extend([name_curr.split(".")[0] for name_curr in curr_salmon["Name"]])
    transcript_counts.extend(curr_salmon["NumReads"])
    sample_names.extend([folder_names[curr_ind].split("_")[-2]] * len(curr_salmon.index))

names_to_reads = pd.DataFrame({"TranscriptNames": transcript_names, "NumReads": transcript_counts, "SampleName": sample_names})
if ".csv" in config["names_to_reads"]:
    names_to_reads.to_csv(path_or_buf = os.path.join(config["reference"],config["names_to_reads"]), sep = "\t")
else:
    names_to_reads.to_csv(path_or_buf = os.path.join(config["reference"],"namestoreads.csv"), sep = "\t")
    config['names_to_reads'] = "namestoreads.csv"

    with open('config.yaml', 'w') as f:
        yaml.dump(config, f)
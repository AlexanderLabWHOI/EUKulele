import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('--out_prefix')
parser.add_argument('--output_dir')
args = parser.parse_args()

with open("config.yaml", 'r') as configfile:
    config = yaml.safe_load(configfile)
return config

results_frame = dict()
met_dir = os.path.join(config["output"],"METs")
samples_dir = os.path.join(config["samples"],"METs")
nucle_extension = config["nucleotide_extension"]
prot_extension = config["protein_extension"]
use_counts = config["use_salmon_counts"]

### READ IN RESULTS FILES FROM MET DIR THAT FIT SAMPLE SPEC FROM CONFIG ###
samples = os.listdir(samples_dir)
for s in samples:
    if ("SH" in s) & ("merged" in s) & ((nucle_extension in s) | (prot_extension in s)):
        file_name = s.split(".")[0] + "-estimated-taxonomy-parallel.out"
        if not os.path.isfile(os.path.join(met_dir,file_name)):
            print("One of the files, " + s + ", in the sample directory did not complete successfully.")
            #sys.exit(1)
        else:
            results_frame[file_name] = pd.read_csv(os.path.join(met_dir,file_name), sep = "\t", index_col=0)

def countClassifs(level, level_hierarchy, name_level, df):
    set_list = set()
    
    classifications = list(df.loc[df["classification_level"] == level]["classification"])
    counts = list(df.loc[df["classification_level"] == level]["counts"])
    match_loc = int(np.where([curr == level for curr in level_hierarchy])[0])
    for curr in range(match_loc + 1,len(level_hierarchy)):
        classification_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["full_classification"])
        correct_index = list(np.where([len(str(cr).split(";")) >= abs(-1-(curr-match_loc)) for cr in classification_curr]))
        
        classification_curr = [classification_curr[cr2] for cr2 in correct_index]
        classifs_curr = [str(cr).split(";")[-1-(curr-match_loc)].strip() for cr in classification_curr]
        counts_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["counts"])
        counts_curr = [classification_curr[cr2] for cr2 in correct_index]
        
        # add to running list
        classifications.extend(classifs_curr)
        counts.extend(counts_curr)
        
    classifications = [cr.strip().strip(""''"]['") for cr in classifications]
    full_list = classifications
    set_list.update(set(classifications))
    
    full_frame = pd.DataFrame({name_level: classifications, "Counts": counts})
    final_frame = full_frame.groupby(name_level)['Counts'].agg(Sum='sum', Count='count')
    return classifications, final_frame
    
def countClassifsNoCounts(level, level_hierarchy, name_level, df):
    set_list = set()
    
    classifications = list(df.loc[df["classification_level"] == level]["classification"])
    match_loc = int(np.where([curr == level for curr in level_hierarchy])[0])
    for curr in range(match_loc + 1,len(level_hierarchy)):
        classification_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["full_classification"])
        classifs_curr = [str(cr).split(";")[-1-(curr-match_loc)].strip() for cr in classification_curr if len(str(cr).split(";")) >= abs(-1-(curr-match_loc))]
        
        # add to running list
        classifications.extend(classifs_curr)
        
    classifications = [cr.strip().strip(""''"]['") for cr in classifications]
    full_list = classifications
    set_list.update(set(classifications))
    final_frame = list(zip(list(set_list), [full_list.count(curr) for curr in list(set_list)]))
    
    return classifications, final_frame

def stripClassifData(df):
    level_hierarchy = ['supergroup','division','class','order','family','genus','species']
    
    return_dict_list = dict()
    return_dict_frame = dict()
    for curr_level in level_hierarchy:
        if use_counts == 1:
            curr_list, curr_counts = countClassifs(curr_level, level_hierarchy, curr_level.capitalize(), df)
            return_dict_list[curr_level] = curr_list
            return_dict_frame[curr_level] = curr_counts
        else:
            curr_list, curr_counts = countClassifsNoCounts(curr_level, level_hierarchy, curr_level.capitalize(), df)
            return_dict_list[curr_level] = curr_list
            return_dict_frame[curr_level] = curr_counts
    return return_dict_list, return_dict_frame


list_results = dict()
frame_results = dict()

# characterizing by major classes 
for curr in results_frame.keys():
    list_results[curr], frame_results[curr] = stripClassifData(results_frame[curr])
    
def makeConcatFrame(curr_df, new_df, level, sample_name):
    if use_counts == 1:
        new_df = pd.DataFrame(new_df.reset_index())
        new_df.columns = [level,"Counts","NumTranscripts"]
    else:
        new_df = pd.DataFrame(new_df)
        new_df.columns = ["Species","NumTranscripts"]
    
    new_df["Sample"] = sample_name
    return pd.concat([curr_df, new_df])

counts_all = dict()
level_hierarchy = ['supergroup','division','class','order','family','genus','species']
for l in level_hierarchy:
    counts_all[l] = pd.DataFrame(columns = ["Species","NumTranscripts","Sample"])
    
for curr in results_frame.keys():
    if "_merged" not in curr:
        continue

    sample_name = curr.split("_")[0]
    
    for l in level_hierarchy:
        curr_df = counts_all[l]
        new_df = frame_results[curr][l]
        counts_all[l] = makeConcatFrame(curr_df, new_df, l.capitalize(), sample_name)
    
for l in counts_all:
    prefix = args.out_prefix
    counts_all[l].to_csv(os.path.join(args.output_dir, prefix + "all_" + l + "_counts.csv"))
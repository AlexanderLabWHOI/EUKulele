import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys
import yaml

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
    if ("merged" in s) & ("estimated" in s) & ((nucle_extension in s) | (prot_extension in s)):
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
    match_loc = np.where(level_hierarchy == level)[0]
    for curr in range(match_loc + 1,len(level_hierarchy)):
        classification_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["full_classification"])
        classifs_curr = [str(cr).split(";")[-1-(curr-match_loc)].strip() for cr in classification_curr]
        counts_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["counts"])
        
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
    match_loc = np.where(level_hierarchy == level)[0]
    for curr in range(match_loc + 1,len(level_hierarchy)):
        classification_curr = list(df.loc[df["classification_level"] == level_hierarchy[curr]]["full_classification"])
        classifs_curr = [str(cr).split(";")[-1-(curr-match_loc)].strip() for cr in classification_curr]
        
        # add to running list
        classifications.extend(classifs_curr)
        
    classifications = [cr.strip().strip(""''"]['") for cr in classifications]
    full_list = classifications
    set_list.update(set(classifications))
    final_frame = list(zip(list(set_list), [full_list.count(curr) for curr in list(set_list)]))
    
    return classifications, final_frame

def stripClassifData(df):
    level_hierarchy = ['supergroup','division','class','order','family','genus','species']
    
    for curr_level in level_hierarchy:
        if use_counts == 1:
            exec(curr_level + "_list," + curr_level + "_counts = " + countClassifs(curr_level, level_hierarchy, curr_level.capitalize(), df))
        else:
            exec(curr_level + "_list," + curr_level + "_counts = " + countClassifsNoCounts(curr_level, level_hierarchy, curr_level.capitalize(), df))
    return species_list, genus_list, family_list, supergroup_counts, species_counts, genus_counts, family_counts


class_frame = dict()
supergroup_counts = dict()
unique_species = dict()
unique_genus = dict()
unique_family = dict()
species_counts = dict()
genus_counts = dict()
family_counts = dict()

# characterizing by major classes 
for curr in results_frame.keys():
    unique_species[curr], unique_genus[curr], unique_family[curr], supergroup_counts[curr], species_counts[curr], genus_counts[curr], family_counts[curr] = stripClassifData(results_frame[curr])
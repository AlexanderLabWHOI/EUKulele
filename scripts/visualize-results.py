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
samples_dir = config["samples"]
nucle_extension = config["nucleotide_extension"]
prot_extension = config["protein_extension"]
use_counts = config["use_salmon_counts"]

### READ IN RESULTS FILES FROM MET DIR THAT FIT SAMPLE SPEC FROM CONFIG ###
samples = os.listdir(samples_dir)
for s in samples:
    if ("merged" in s) & ("estimated" in s) & ((nucle_extension in s) | (prot_extension in s)):
        if !os.path.isfile(s):
            print("Not all of the files in the sample directory completed successfully.")
            sys.exit(1)
        results_frame[s] = pd.read_csv(os.path.join(met_dir,s), sep = "\t", index_col=0)
  

class_frame = dict()
supergroup_counts = dict()
unique_species = dict()
unique_genus = dict()
unique_family = dict()
species_counts = dict()
genus_counts = dict()
family_counts = dict()

def countClassifs(level, level_hierarchy):
    classifications = list(df.loc[df["classification_level"] == level]["classification"])
    full_classification = list(df.loc[df["classification_level"] == level]["full_classification"])
    match_loc = np.where(level_hierarchy == level)[0]
    for curr in range(match_loc + 1,len(level_hierarchy):
        classification 

def stripClassifData(df):
    species = set()
    families = set()
    genuses = set()
    supergroups = set()
    level_hierarchy = ['supergroup','division','class','order','family','genus','species']
    full_class_spec = list(df.loc[df["classification_level"] == "species"]["full_classification"])
    full_class_spec_ct = list(df.loc[df["classification_level"] == "species"]["counts"])
    full_class_genus = list(df.loc[df["classification_level"] == "genus"]["full_classification"])
    full_class_genus_ct = list(df.loc[df["classification_level"] == "genus"]["counts"])
    full_class_family = list(df.loc[df["classification_level"] == "family"]["full_classification"])
    full_class_family_ct = list(df.loc[df["classification_level"] == "family"]["counts"])
    
    full_class_spec = [cr.strip().strip(""''"]['") for cr in full_class_spec]
    full_class_genus = [cr.strip().strip(""''"]['") for cr in full_class_genus]
    full_class_family = [cr.strip().strip(""''"]['") for cr in full_class_family]
    
    specie = list(df.loc[df["classification_level"] == "species"]["classification"])
    specie_ct = full_class_spec_ct
    genus = [str(cr).split(";")[len(str(cr).split(";"))-2].strip() for cr in full_class_spec] + \
            list(df.loc[df["classification_level"] == "genus"]["classification"])
    genus_ct = full_class_spec_ct + full_class_genus_ct # if len(cr.split(";")) >= 3] + \
    family = [str(cr).split(";")[len(cr.split(";"))-3].strip() for cr in full_class_spec] + \
            [str(cr).split(";")[len(str(cr).split(";"))-2].strip() for cr in full_class_genus] + \
            list(df.loc[df["classification_level"] == "family"]["classification"])
    family_ct = full_class_spec_ct + full_class_genus_ct + full_class_family_ct
    
    supergroup = [str(cr).split(";")[0].strip() for cr in list(df["full_classification"])] #+ \
            #list(df.loc[df["classification_level"] == "supergroup"]["classification"])
    
    # Add to supergroup set
    full_list = supergroup
    supergroups.update(set(supergroup))
    full_frame = pd.DataFrame({"Supergroup": supergroup, "Counts": list(df["counts"])})
    #supergroup_counts = list(zip(list(supergroups), [full_list.count(curr) for curr in list(supergroups)]))
    supergroup_counts = full_frame.groupby('Supergroup')['Counts'].agg(Sum='sum', Count='count')
    
    # Add to families set
    full_list = family
    families.update(set(family))
    full_frame = pd.DataFrame({"Family": family, "Counts": list(family_ct)})
    #families_counts = list(zip(list(families), [full_list.count(curr) for curr in list(families)]))
    families_counts = full_frame.groupby('Family')['Counts'].agg(Sum='sum', Count='count')
    
    # Add to genus set
    full_list = genus
    genuses.update(set(genus))
    full_frame = pd.DataFrame({"Genus": genus, "Counts": list(genus_ct)})
    #genuses_counts = list(zip(list(genuses), [full_list.count(curr) for curr in list(genuses)]))
    genuses_counts = full_frame.groupby('Genus')['Counts'].agg(Sum='sum', Count='count')
    
    # Add to species set
    full_list = specie
    species.update(set(specie))
    full_frame = pd.DataFrame({"Species": specie, "Counts": list(specie_ct)})
    #species_counts = list(zip(list(species), [full_list.count(curr) for curr in list(species)]))
    species_counts = full_frame.groupby('Species')['Counts'].agg(Sum='sum', Count='count')
    
    return species, genuses, families, supergroup_counts, species_counts, genuses_counts, families_counts
    
# characterizing by major classes 
for curr in results_frame.keys():
    unique_species[curr], unique_genus[curr], unique_family[curr], supergroup_counts[curr], species_counts[curr], genus_counts[curr], family_counts[curr] = stripSpeciesGenusAndOrder(results_frame[curr])
    
    class_frame[curr] = sum_results
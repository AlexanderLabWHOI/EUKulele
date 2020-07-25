import os
import sys
import subprocess

def readBuscoFile(individual_or_summary, busco_file, organisms, organisms_taxonomy):
    if individual_or_summary == "individual":
        if (busco_file != "") & (os.path.isfile(busco_file)):
            busco_file_read = read.csv(busco_file, sep = "\t")
            organisms = list(busco_file_read.iloc[:,0])
            organisms_taxonomy = list(busco_file_read.iloc[:,1])
            print("Organisms and their taxonomy levels for BUSCO analysis were read from file.")
            
            if (len(organisms) != len(organisms_taxonomy)):
                print("Organisms and taxonomic specifications for BUSCO analysis do not contain the same number of entries. " + \
                      "Please revise such that each organism flagged for BUSCO analysis also includes its original taxonomic " + \
                      "level.")
                sys.exit(1)
        else:
            print("No BUSCO file specified/found; using argument-specified organisms and taxonomy for BUSCO analysis.")
            
    return organisms, organisms_taxonomy
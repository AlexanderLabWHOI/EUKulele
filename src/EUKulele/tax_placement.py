'''
Script for handling taxonomic estimation.
'''

#!/usr/bin/env python
import os
import sys
import chardet
import argparse
import multiprocessing
from joblib import Parallel, delayed
import ujson
import yaml
import pandas as pd
import numpy as np

import EUKulele

def tax_placement(pident, tax_cutoffs):
    ''' Decide which level of taxonomy is appropriate. '''

    if pident >= tax_cutoffs['species']:
        out = 'species'
        level = 8
    elif pident >= tax_cutoffs['genus']:
        out = 'genus'
        level = 7
    elif pident >= tax_cutoffs['family']:
        out = 'family'
        level = 6
    elif pident >= tax_cutoffs['order']:
        out = 'order'
        level = 5
    elif pident >= tax_cutoffs['class']:
        out = 'class'
        level = 4
    elif pident >= tax_cutoffs['division']:
        out = 'division'
        level = 3
    elif pident >= tax_cutoffs['supergroup']:
        out = 'supergroup'
        level = 2
    elif pident >= tax_cutoffs['domain']:
        out = 'domain'
        level = 1
    else:
        out = 'unclassified'
        level = 0
        
    return out, level

def read_in_taxonomy(infile):
    ''' Load taxonomy table. '''

    with open(infile, 'rb') as f:
        result = chardet.detect(f.read())
    tax_out = pd.read_csv(infile, sep='\t',encoding=result['encoding'])
    classes = ['domain','supergroup','division','class','order',
               'family','genus','species']
    classes_out = []
    for c in tax_out.columns:
        if c.lower() in classes:
            if (np.issubdtype(tax_out[str(c)].dtypes, np.number)) | \
               (np.issubdtype(tax_out[str(c)].dtypes, np.float_)):
                tax_out = tax_out.loc[:,~(tax_out.columns == c)]
            classes_out.append(c.lower())
    tax_out.columns = tax_out.columns.str.lower()
    tax_out = tax_out.set_index('source_id')
    return tax_out, classes_out

def read_in_tax_cutoffs(yamlfile):
    ''' Read taxonomic cutoff YAML file.'''

    with open(yamlfile, 'r') as stream:
        co_out = yaml.safe_load(stream)
    return co_out

def read_in_protein_map(protjson):
    ''' Read protein map JSON file. '''

    with open(protjson, 'rb') as f:
        pout = ujson.load(f)
    return pout

def gen_dict(tax_table,classes):
    ''' Create classication dictionary for easy lookup. '''

    tax_table["Classification"] = ""
    tax_table = tax_table.loc[:,~tax_table.columns.duplicated()]
    for c in classes:
        tax_table.loc[tax_table[c] != tax_table[c], c] = "No " + c
        if all([str(curr).lower() != "nan" for curr in list(tax_table[c])]):
            if str(tax_table["Classification"][0]) != "":
                tax_table["Classification"] = tax_table["Classification"] +\
                                              ";" + tax_table[c]
            else:
                tax_table["Classification"] = tax_table["Classification"] +\
                                              tax_table[c]
    return dict(zip(tax_table.index, tax_table["Classification"]))

def gen_reads_dict(names_to_reads):
    names_to_reads = pd.read_csv(names_to_reads,header=0,sep="\t")
    return dict(zip(names_to_reads["TranscriptNames"],
                    names_to_reads["NumReads"]))

def lca(full_classifications,classes):
    ''' Last Common Ancestor (LCA) algorithm. '''

    full_classifications_split = [[str(subtax).strip() for
                                   subtax in curr.split(";")] for
                                  curr in full_classifications]
    length_classes = [len(curr) for curr in full_classifications_split]
    if len(set(length_classes)) != 1:
        print("Error: not all classifications at at the same taxonomic level.",
              flush = True)
        sys.exit(1)
    for l in reversed(range(length_classes[0])):
        set_classifications = [curr[l] for curr in full_classifications_split]
        if len(set(set_classifications)) == 1:
            return classes[l], set_classifications[0], \
                   "; ".join(full_classifications_split[0][0:(l+1)])
    return "","","" # if there are no common ancestors

def match_maker(dd, consensus_cutoff, tax_dict, use_counts, tax_cutoffs, classes):
    ''' Manages decision between multiple matches. '''

    ambiguous = 0 # we assume unambiguous
    md = dd.pident.max()
    transcript_name = set(list(dd["qseqid"]))
    if len(transcript_name) > 1:
        print("More than 1 transcript name included in the group.", flush = True)
    transcript_name = list(transcript_name)[0]
    ds = list(set(dd[dd.pident==md]['ssqid_TAXID']))
    counts = list(set(dd[dd.pident==md]['counts']))
    if len(counts) >= 1:
        chosen_count = counts[0]
    else:
        chosen_count = 0
    assignment, level = tax_placement(md, tax_cutoffs)
                        # most specific taxonomic level assigned
    if len(ds)==1:
        if ds[0] not in tax_dict:
            return pd.DataFrame(columns=['transcript_name','classification_level',
                                         'full_classification','classification',
                                         'max_pid','counts','ambiguous'])
        full_classification = str(tax_dict[ds[0]]).split(";")[0:level]
        best_classification = full_classification[len(full_classification) - 1]
                            # the most specific taxonomic level we can classify by
        full_classification = '; '.join(full_classification)
                              # the actual assignments based on that
    else:
        classification_0 = []
        full_classification_0 = []
        for d in ds:
            if d not in tax_dict:
                return(pd.DataFrame(columns=['transcript_name','classification_level',
                                             'full_classification','classification',
                                             'max_pid','counts','ambiguous']))
            d_full_class = str(tax_dict[d]).split(";")[0:level]
            classification_0.append(d_full_class[len(d_full_class) - 1])
                        # the most specific taxonomic level we can classify by
            full_classification_0.append('; '.join(d_full_class))
                        # the actual assignments based on that
        entries = list(set(full_classification_0))
        if len(entries) == 1:
            best_classification = classification_0[0]
            full_classification = full_classification_0[0]
        else:
            ambiguous = 1
            best_frac = 0
            best_one_class = 0
            best_full_class = 0
            for e in entries:
                curr_frac = len(np.where(full_classification_0 == e)) / \
                            len(full_classification_0)
                if (isinstance(curr_frac, float)) & (curr_frac > best_frac):
                    best_frac = curr_frac
                    best_full_class = e
                    best_one_class = str(e.split(";")[len(e.split(";")) - 1]).strip()
            if best_frac >= consensus_cutoff:
                best_classification = best_one_class
                full_classification = best_full_class
            else:
                assignment, best_classification, \
                        full_classification = lca(full_classification_0,classes)

    if use_counts == 1:
        return pd.DataFrame([[transcript_name, assignment,\
                              full_classification, best_classification, md,\
                              chosen_count, ambiguous]],\
                       columns=['transcript_name','classification_level', 'full_classification',
                                'classification', 'max_pid', 'counts', 'ambiguous'])
    else:
        return pd.DataFrame([[transcript_name, assignment, full_classification,\
                              best_classification, md, ambiguous]],\
                       columns=['transcript_name', 'classification_level', 'full_classification',
                                'classification', 'max_pid', 'ambiguous'])

def apply_parallel(grouped_data, match_maker, consensus_cutoff, tax_dict, use_counts, tax_cutoffs, classes):
    resultdf = Parallel(n_jobs=multiprocessing.cpu_count(),
                        prefer="threads")(delayed(match_maker)(group,
                                                               consensus_cutoff,
                                                               tax_dict,
                                                               use_counts,
                                                               tax_cutoffs,
                                                               classes)
                                          for name, group in grouped_data)
    return pd.concat(resultdf)

def classify_taxonomy_parallel(df, tax_dict, namestoreads, pdict,
                               consensus_cutoff, tax_cutoffs, classes):
    ''' Parallel implementation of the taxonomic classication process. '''

    chunksize = 2 * 10 ** 6
    counter = 0

    ## Return an empty dataframe if no matches made ##
    if os.stat(str(df)).st_size == 0:
        if namestoreads != 0:
            return pd.DataFrame(columns=['transcript_name','classification_level',
                                         'full_classification',
                                'classification', 'max_pid', 'counts', 'ambiguous'])
        else:
            return pd.DataFrame(columns=['transcript_name','classification_level',\
                                         'full_classification',
                                'classification', 'max_pid', 'ambiguous'])

    for chunk in pd.read_csv(str(df), sep = '\t', header = None, chunksize=chunksize):
        chunk.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                         'gapopen', 'qstart','qend', 'sstart', 'send',
                         'evalue', 'bitscore']
        chunk.sseqid = [curr.replace(".","N") for curr in chunk.sseqid]
        chunk['ssqid_TAXID']=chunk.sseqid.map(pdict)
        chunk = chunk[['qseqid','pident', 'evalue', 'bitscore', 'ssqid_TAXID']]
        if namestoreads != 0:
            chunk['counts']=[namestoreads[curr.split(".")[0].split(":")[0]] \
                             if curr.split(".")[0].split(":")[0] in \
                             namestoreads else 0 for curr in chunk.qseqid]
            use_counts = 1
        else:
            chunk['counts'] = [0] * len(chunk.qseqid)
            # if no reads dict, each count is just assumed to be 0 and isn't recorded later
            use_counts = 0

        if counter == 0:
            outdf = apply_parallel(chunk.groupby('qseqid'), match_maker,
                                   consensus_cutoff, tax_dict, use_counts, tax_cutoffs, classes)
        else:
            outdf = pd.concat([outdf, apply_parallel(chunk.groupby('qseqid'),
                                                     match_maker, consensus_cutoff, tax_dict,
                                                     use_counts, tax_cutoffs, classes)], axis = 0)
        counter = counter + 1
    return outdf

def place_taxonomy(tax_file,cutoff_file,consensus_cutoff,prot_map_file,
                   use_counts,names_to_reads,diamond_file,outfile,rerun):
    ''' Find predicted taxonomy using alignment matches. '''

    if (os.path.isfile(outfile)) & (not rerun):
        print("Taxonomic placement already complete at", outfile + "; will not re-run step.")
        return pd.read_csv(outfile, sep = "\t")
 
    tax_table, classes = read_in_taxonomy(tax_file)
    tax_cutoffs = read_in_tax_cutoffs(os.path.join(os.path.dirname(\
        os.path.realpath(__file__)), "static", cutoff_file))
    pdict = read_in_protein_map(prot_map_file)
    tax_dict = gen_dict(tax_table,classes)
    consensus_cutoff = float(consensus_cutoff)
    if int(use_counts) == 1:
        reads_dict = gen_reads_dict(names_to_reads)
        classification_df = classify_taxonomy_parallel(diamond_file, tax_dict, reads_dict,
                                                       pdict, consensus_cutoff, tax_cutoffs, classes)
    else:
        classification_df = classify_taxonomy_parallel(diamond_file, tax_dict, 0, pdict,
                                                       consensus_cutoff, tax_cutoffs, classes)
    classification_df.to_csv(outfile, sep='\t')
    return outfile

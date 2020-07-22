#!/usr/bin/env python
import ujson
import pandas as pd
import yaml
import chardet
import argparse
import numpy as np
import os

level_dict = {'class':['supergroup','division','class'], 'order':['supergroup','division','class', 'order'], 'family':['supergroup','division','class', 'order', 'family'], 'genus': ['supergroup','division','class', 'order', 'family','genus'], 'species':['supergroup','division','class', 'order', 'family','genus', 'species']}
levels = ['supergroup','division','class','order','family','genus','species']

# bug with splitting based on the new taxonomy estimation
def split_taxonomy(protein_estimate):
    outdf = protein_estimate.copy()
    split_tax_nan = pd.DataFrame(protein_estimate.full_classification.str.split('; '))
    new_col_dict = dict()
    for l in levels:
        new_col_dict[str(l)] = [""] * len(outdf.index)
    for cl in range(len(split_tax_nan.index)):
        lineage = list(split_tax_nan.loc[:,'full_classification'])[cl]
        for i,l in enumerate(levels):
            if isinstance(lineage,list):
                if (len(lineage)>i):
                    new_col_dict[l][cl] = lineage[i]
                else:
                    new_col_dict[l][cl] = np.nan
            else:
                new_col_dict[l][cl] = np.nan
    for l in levels:
        outdf[str(l)] = new_col_dict[l]
    return(outdf)

def create_tax_dictionary(split_taxonomy_df):
    tax_dict = {}
    total_len=len(split_taxonomy_df)
    for l in levels:
        tax_dict[l]={}
        column_sum = split_taxonomy_df.groupby(l)[l].sum()
        column_count = split_taxonomy_df.groupby(l)[l].count()
        norm_count=column_count/total_len
        tax_dict[l]=norm_count
    return(tax_dict)

def get_max_levels(tax_dict):
    max_df = pd.DataFrame(index=levels, columns=['max_taxa','percent_id'])
    for key in tax_dict:
        highest_tax = tax_dict[key][tax_dict[key]==tax_dict[key].max()].index
        max_val = tax_dict[key].max()
        if len(highest_tax)>0:
            max_df.loc[key]=highest_tax[0], max_val
        else:
            max_df.loc[key]=np.nan, max_val
    return(max_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--estimated-taxonomy-file')
    parser.add_argument('--out-prefix')
    parser.add_argument('--outdir')
    parser.add_argument('--max-out-dir')
    args = parser.parse_args()
    estimated_tax = pd.read_csv(args.estimated_taxonomy_file, sep='\t', index_col=0)
    split_taxonomy_df = split_taxonomy(estimated_tax)
    tax_dict = create_tax_dictionary(split_taxonomy_df)
    max_df = get_max_levels(tax_dict)
    if not os.path.exists(args.outdir):
        try:
            os.mkdir(args.outdir)
        except:
            pass
    if not os.path.exists(args.max_out_dir):
        try:
            os.mkdir(args.max_out_dir)
        except:
            pass
    for l in levels:
        tax_dict[l].to_csv(os.path.join(args.outdir, args.out_prefix+'.'+l), header=False, sep='\t')
    max_df.to_csv(os.path.join(args.max_out_dir, args.out_prefix + '-max-level.csv'), sep='\t')

'''
Function file containing information on how to calculate
maximum taxonomic level among MAGs.
'''

import os
import argparse
import pandas as pd
import numpy as np

level_dict = {'class':['domain','supergroup','division','class'],
              'order':['domain','supergroup','division','class','order'],
              'family':['domain','supergroup','division','class','order','family'],
              'genus': ['domain','supergroup','division','class','order','family','genus'],
              'species':['domain','supergroup','division','class','order','family','genus','species'],}
levels = ['domain','supergroup','division','class','order','family','genus','species']

# bug with splitting based on the new taxonomy estimation
def split_taxonomy(protein_estimate, level_hierarchy):
    ''' Split the taxonomy based on the matched protein.'''
    outdf = protein_estimate.copy()
    split_tax_nan = pd.DataFrame(protein_estimate.full_classification.str.split(';'))
    for col in split_tax_nan.columns:
        if split_tax_nan[col].dtype == "object":
            split_tax_nan[col] = [[curr2.strip() for curr2 in curr] for curr in split_tax_nan[col]]
            
    new_col_dict = dict()
    for curr_level in level_hierarchy:
        new_col_dict[str(curr_level)] = [""] * len(outdf.index)
    for curr_cl in range(len(split_tax_nan.index)):
        lineage = list(split_tax_nan.loc[:,'full_classification'])[curr_cl]
        for i,curr_level in enumerate(level_hierarchy):
            if curr_level in new_col_dict:
                if isinstance(lineage,list):
                    if len(lineage) > i:
                        new_col_dict[curr_level][curr_cl] = lineage[i]
                    else:
                        new_col_dict[curr_level][curr_cl] = np.nan
                else:
                    new_col_dict[curr_level][curr_cl] = np.nan
    for curr_level in level_hierarchy:
        if curr_level in new_col_dict:
            outdf[str(curr_level)] = new_col_dict[curr_level]
    return outdf

def create_tax_dictionary(split_taxonomy_df,level_hierarchy):
    ''' Create a Python dictionary to quickly check taxonomy.'''
    tax_dict = {}
    total_len=len(split_taxonomy_df)
    for curr_level in level_hierarchy:
        if curr_level in split_taxonomy_df.columns:
            tax_dict[curr_level]={}
            # column_sum = split_taxonomy_df.groupby(l)[curr_level].sum()
            column_count = split_taxonomy_df.groupby(curr_level)[curr_level].count()
            norm_count=column_count/total_len
            tax_dict[curr_level]=norm_count
    return tax_dict

def get_max_levels(tax_dict,level_hierarchy):
    '''Get the most confident taxonomic assignment.'''
    max_df = pd.DataFrame(index=level_hierarchy, columns=['max_taxa','percent_id'])
    for key in tax_dict:
        highest_tax = tax_dict[key][tax_dict[key]==
                                    tax_dict[key].max()].index
        max_val = tax_dict[key].max()
        if len(highest_tax)>0:
            max_df.loc[key]=highest_tax[0], max_val
        else:
            max_df.loc[key]=np.nan, max_val
    return max_df

def magStats(args=None):
    '''Create dataframe summarizing output.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('--estimated-taxonomy-file')
    parser.add_argument('--out-prefix')
    parser.add_argument('--outdir')
    parser.add_argument('--max-out-dir')
    parser.add_argument('--level-hierarchy')
    
    if args is not None:
        args = parser.parse_args(args)
    else:
        args = parser.parse_args()
        
    level_hierarchy=args.level_hierarchy.split(" ")
    os.system("mkdir -p " + args.max_out_dir)
    os.system("mkdir -p " + args.outdir)
    estimated_tax = pd.read_csv(args.estimated_taxonomy_file, sep='\t', index_col=0)
    split_taxonomy_df = split_taxonomy(estimated_tax,level_hierarchy)
    #split_taxonomy_df[level_hierarchy] = split_taxonomy_df.str.split(full_classification,sep=";")
    tax_dict = create_tax_dictionary(split_taxonomy_df,level_hierarchy)
    max_df = get_max_levels(tax_dict,level_hierarchy)
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
    for curr_level in args.level_hierarchy:
        if curr_level in tax_dict:
            tax_dict[curr_level].to_csv(os.path.join(args.outdir,
                                            args.out_prefix + '.' +\
                                                     curr_level),
                               header=False, sep='\t')
    max_df.to_csv(os.path.join(args.max_out_dir,
                               args.out_prefix +\
                               '-max-level.csv'), sep='\t')
    return 0

if __name__ == "__main__":
    magStats()

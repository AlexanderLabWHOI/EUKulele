#!/usr/bin/env python
import ujson
import pandas as pd
import yaml
import chardet
import argparse
import multiprocessing
from joblib import Parallel, delayed

def tax_placement(pident):
    if pident >= tax_cutoffs['species']:
        out = 'species'; level = 6;
    elif pident >= tax_cutoffs['genus']:
        out = 'genus'; level = 5;
    elif pident >= tax_cutoffs['family']:
        out = 'family'; level = 4;
    elif pident >= tax_cutoffs['order']:
        out = 'order'; level = 3;
    elif pident < tax_cutoffs['order']:
        out = 'class'; level = 2;
    return out, level

def read_in_taxonomy(infile):
    with open(infile, 'rb') as f:
        result = chardet.detect(f.read())
    tax_out = pd.read_csv(infile, sep='\t',encoding=result['encoding'])
    tax_out.columns = tax_out.columns.str.lower()
    tax_out = tax_out.set_index('source_id')
    return tax_out

def read_in_tax_cutoffs(yamlfile):
    with open(yamlfile, 'r') as stream:
        co_out = yaml.safe_load(stream)
    return co_out

def read_in_protein_map(protjson):
    with open(protjson, 'rb') as f:
        pout = ujson.load(f)
    return pout

def read_in_diamond_file(dfile, pdict):
    dfout =  pd.read_csv(dfile, sep = '\t', header = None)
    dfout.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    dfout['ssqid_TAXID']=dfout.sseqid.map(pdict)
    return dfout

def gen_dict(tax_table):
    classes = ['supergroup','division','class','order','family','genus','species']
    tax_table["Classification"] = ""
    for c in classes:
        if str(tax_table["Classification"][0]) != "":
            tax_table["Classification"] = tax_table["Classification"] + ";" + tax_table[c]
        else:
            tax_table["Classification"] = tax_table["Classification"] + tax_table[c]
    return(dict(zip(tax_table.index, tax_table["Classification"])))

def match_maker(dd):
    ambiguous = 0 # we assume unambiguous
    md = dd.pident.max()
    ds = list(set(dd[dd.pident==md]['ssqid_TAXID']))
    assignment, level = tax_placement(md) # most specific taxonomic level assigned
    if len(ds)==1:
        full_classification = tax_dict[ds[0]].split(";")[0:level]
        best_classification = full_classification[len(full_classification) - 1] # the most specific taxonomic level we can classify by
        full_classification = '; '.join(full_classification) # the actual assignments based on that
    else:
        classification_0 = set()
        full_classification_0 = set()
        for d in ds:
            d_full_class = tax_dict[d].split(";")[0:level]
            classification_0.add(d_full_class[len(d_full_class) - 1]) # the most specific taxonomic level we can classify by
            full_classification_0.add('; '.join(d_full_class)) # the actual assignments based on that
            if len(classification_0) == 1:
                best_classification = list(set(classification_0))[0]
                full_classification = list(set(full_classification_0))[0]
            else:
                ambiguous = 1
                break # we just end up picking the first added if we have multiple classifications.
                # we probably want to implement this differently in the future. For now we mark as "ambiguous"
    return pd.DataFrame([[assignment, full_classification, best_classification, md, ambiguous]],\
                   columns=['classification_level', 'full_classification', 'classification', 'max_pid', 'ambiguous'])

def apply_parallel(grouped_data, match_maker):
    resultdf = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(match_maker)(group) for name, group in grouped_data)
    return pd.concat(resultdf)

def classify_taxonomy_parallel(df, tax_dict):
    outdf = apply_parallel(df.groupby('qseqid'), match_maker)
    return outdf

# create a dictionary for all of the mmetsp and then just split by ";" and then take the top X based on the tax class level.
# time the difference between the new function and the original.
def classify_taxonomy(df, tax_dict):
    level_dict = {'class':['supergroup','division','class'],
                       'order':['supergroup','division','class', 'order'],
                       'family':['supergroup','division','class', 'order', 'family'],
                       'genus': ['supergroup','division','class', 'order', 'family','genus'],
                        'species':['supergroup','division','class', 'order', 'family','genus', 'species']}

    outdf = pd.DataFrame(columns = ['classification_level', 'full_classification', 'classification', 'max_pid', 'ambiguous'])
    for t,dd in df.groupby('qseqid'):
        ambiguous = 0 # we assume unambiguous
        md = dd.pident.max()
        ds = list(set(dd[dd.pident==md]['ssqid_TAXID']))
        assignment, level = tax_placement(md) # most specific taxonomic level assigned
        if len(ds)==1:
            full_classification = tax_dict[ds[0]].split(";")[0:level]
            best_classification = full_classification[len(full_classification) - 1] # the most specific taxonomic level we can classify by
            full_classification = '; '.join(full_classification) # the actual assignments based on that
        else:
            classification_0 = set()
            full_classification_0 = set()
            for d in ds:
                d_full_class = tax_dict[d].split(";")[0:level]
                classification_0.add(d_full_class[len(d_full_class) - 1]) # the most specific taxonomic level we can classify by
                full_classification_0.add('; '.join(d_full_class)) # the actual assignments based on that
                if len(classification_0) == 1:
                    best_classification = list(classification_0)[0]
                    full_classification = list(full_classification_0)[0]
                else:
                    ambiguous = 1
                    break # we just end up picking the first added if we have multiple classifications.
                    # we probably want to implement this differently in the future. For now we mark as "ambiguous"
        outdf.loc[t] =  [assignment, full_classification, best_classification, md, ambiguous]
    return outdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tax_file')
    parser.add_argument('--cutoff_file')
    parser.add_argument('--prot_map_file')
    parser.add_argument('--diamond_file')
    parser.add_argument('--outfile')
    parser.add_argument('--method')
    args = parser.parse_args()
    tax_table = read_in_taxonomy(args.tax_file)
    tax_cutoffs = read_in_tax_cutoffs(args.cutoff_file)
    pdict = read_in_protein_map(args.prot_map_file)
    diamond_df = read_in_diamond_file(args.diamond_file, pdict)
    tax_dict = gen_dict(tax_table)
    if args.method == "parallel":
        classification_df = classify_taxonomy_parallel(diamond_df, tax_dict)
    else:
        classification_df = classify_taxonomy(diamond_df, tax_dict)
    classification_df.to_csv(args.outfile, sep='\t')

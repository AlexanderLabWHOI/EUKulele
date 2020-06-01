#!/usr/bin/env python
import ujson
import pandas as pd
import yaml
import chardet
import argparse



def tax_placement(pident):
    if pident >= tax_cutoffs['species']:
        out = 'species'
    elif pident >= tax_cutoffs['genus']:
        out = 'genus'
    elif pident >= tax_cutoffs['family']:
        out = 'family'
    elif pident >= tax_cutoffs['order']:
        out = 'order'
    elif pident < tax_cutoffs['order']:
        out = 'class'
    return out

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

def classify_taxonomy(df, tax_table):
    level_dict = {'class':['supergroup','division','class'],
                       'order':['supergroup','division','class', 'order'],
                       'family':['supergroup','division','class', 'order', 'family'],
                       'genus': ['supergroup','division','class', 'order', 'family','genus'],
                        'species':['supergroup','division','class', 'order', 'family','genus', 'species']}

    outdf = pd.DataFrame(columns = ['classification_level', 'full_classification', 'classification', 'max_pid'])
    for t,dd in df.groupby('qseqid'):
        md = dd.pident.max()
        ds = list(set(dd[dd.pident==md]['ssqid_TAXID']))
        assignment = tax_placement(md)
        full_assignment = level_dict[assignment]
        if len(ds)==1:
            classification = tax_table.loc[ds[0], assignment]
            full_classification = tax_table.loc[ds[0], full_assignment]
            full_classification='; '.join(list(full_classification))
        else:
            classification_0 = []
            full_classification_0=[]
            for d in ds:
                classification_0.append(tax_table.loc[ds[0], assignment])
                fc = tax_table.loc[ds[0], full_assignment]
                full_classification_0.append('; '.join(list(fc)))
                if len(set(classification_0)) ==1:
                    classification = list(set(classification_0))[0]
                    full_classification=list(set(full_classification_0))[0]
                else:
                    pass
        outdf.loc[t] =  [assignment, full_classification, classification, md]
    return outdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tax_file')
    parser.add_argument('--cutoff_file')
    parser.add_argument('--prot_map_file')
    parser.add_argument('--diamond_file')
    parser.add_argument('--outfile')
    args = parser.parse_args()
    tax_table = read_in_taxonomy(args.tax_file)
    tax_cutoffs = read_in_tax_cutoffs(args.cutoff_file)
    pdict = read_in_protein_map(args.prot_map_file)
    diamond_df = read_in_diamond_file(args.diamond_file, pdict)
    classification_df = classify_taxonomy(diamond_df, tax_table)
    classification_df.to_csv(args.outfile, sep='\t')

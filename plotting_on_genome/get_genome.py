#!/usr/bin/env python

# Adapted from loadGenomeInfoFromNCBI.py

from Bio import Entrez
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm


def get_all_contigs(genome_ids, fname_prefix):
    records = [Entrez.efetch(db="nucleotide", id=x, 
                             rettype="fasta", retmode="text").read() 
               for x in tqdm(genome_ids)]

    filename = f"{fname_prefix}_combined_contigs.fasta"
    with open(filename, "w") as fh:
        fh.writelines(records)

def get_all_feature_tables(genome_ids, fname_prefix):
    for genome_id in tqdm(genome_ids):
        record = Entrez.efetch(db="nucleotide", 
                               id=genome_id, 
                               rettype="ft", 
                               retmode="text").read()[:-1]
        
        filename = f"{fname_prefix}_{genome_id}.ft"
        with open(filename, 'w') as f:
            f.write(record)

def open_feature_table(filename):
    data = []
    
    with open(filename) as fh:
        row = {}
        cur_accession = fh.readline().split("|")[1]
        for i, line in enumerate(fh):
            if not line.startswith("\t"):
                if row: data.append(row)
                line = line.split("\t")
                row = {"genomic_accession": cur_accession, 
                       "start": int(line[0].strip("<>")), 
                       "end": int(line[1].strip("<>")), 
                       "type": line[2].strip("\n") if len(line) > 2 else None}
            else:
                key, val = line.split("\t")[-2:]
                row.update({key: val.strip("\n")})
    
    df_res =  pd.DataFrame(data)
    
    cols = ["genomic_accession", "start", "end", "locus_tag", "old_locus_tag"]
    df_res = df_res.reindex(df_res.columns.union(cols, sort=False), 
                            axis=1, fill_value=np.nan)
    
    # Add new column to indicate strand direction
    cond = (df_res["end"] > df_res["start"])
    df_res["strand"] = cond.apply(lambda x : "+" if x else "-")
    
    # If there is no old tag, use new
    df_res["old_locus_tag"].fillna(df_res["locus_tag"], inplace=True)
    return df_res

def main(search_term, email, fname_prefix):
    Entrez.email = email

    handle = Entrez.esearch(db='nucleotide', term=search_term, retmax='100')
    res = Entrez.read(handle)

    genome_ids = res['IdList']

    get_all_contigs(genome_ids, fname_prefix)
    get_all_feature_tables(genome_ids, fname_prefix)

    df_annot = pd.concat([open_feature_table(f"{fname_prefix}_{x}.ft") 
                          for x in genome_ids], 
                         ignore_index=True)
    df_annot.to_csv(f"{fname_prefix}_combined_locus_tags.csv", index=False)



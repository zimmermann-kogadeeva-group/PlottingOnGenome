#!/usr/bin/env python

from Bio import Entrez
import numpy as np
import pandas as pd


def get_row(name, section):
    
    old_cols = ["s.start", "s.end", "subject"]
    row = {"Name": name, 
           "f_start": -1, "f_end": -1, "f_alignment_subject": None,
           "r_start": -1, "r_end": -1, "r_alignment_subject": None}
    
    f_cond = section["query"].str.endswith("_F")
    if f_cond.any():
        new_cols = {"s.start": "f_start", "s.end": "f_end", 
                    "subject": "f_alignment_subject"}
        # TODO: how to pick which row is used - first row is used atm
        row.update(section[f_cond].iloc[0][old_cols]
                   .rename(index=new_cols)
                   .to_dict())
    
    r_cond = section["query"].str.endswith("_R")
    if r_cond.any():
        new_cols = {"s.start": "r_start", "s.end": "r_end", 
                    "subject": "r_alignment_subject"}
        row.update(section[r_cond].iloc[0][old_cols]
                   .rename(index=new_cols)
                   .to_dict())        
    
    return row

def get_percent(row, start, end):
    return (100 * 
            (min(row['end'], end) - max(row['start'], start)) / 
            (row['end'] - row['start']))

def get_gene_text(df_annot, start, end, alignment_subject):
    cond1 = (start < df_annot["start"]) & (df_annot["start"] < end)
    cond2 = (start < df_annot["end"]) & (df_annot["end"] < end)
    cond3 = (df_annot["genomic_accession"] == alignment_subject)
    subset = df_annot[(cond1 | cond2) & cond3].drop_duplicates(subset="locus_tag")

    percents = subset.apply(lambda x : get_percent(x, start, end), axis=1)
    gene_text = [f"{x} ({y}, {z:.2f})" for x, y, z in 
                 zip(subset["locus_tag"], subset["strand"], percents)]
    return ", ".join(gene_text)

def main(locus_file, blast_file, outputfile):

    df_annot = pd.read_csv(locus_file)

    cols = ['query', 'subject', 'identity', 'alignment_length', 'mismatches', 
            'gap_opens', 'q.start', 'q.end', 's.start','s.end', 'evalue', 
            'bit score']
    alignment = pd.read_csv(blast_file, sep = '\t', names=cols)
    alignment["query"] = alignment["query"].str.replace("_and_pZE21_R", "")

    alignment["direction"] = (alignment["query"]
                              .apply(lambda x : "F" if x.endswith("_F") else "R"))
    alignment["name"] = (alignment["query"]
                         .str.replace("_F", "")
                         .str.replace("_R", ""))

    res = pd.DataFrame([get_row(name, x) for name, x 
                        in alignment.groupby(by="name")])
    res["r_alignment_subject"].fillna(res["f_alignment_subject"], inplace=True)
    res["f_alignment_subject"].fillna(res["r_alignment_subject"], inplace=True)
    
    mock_insert_length = 4000
    
    cols = ["f_start", "f_end"]
    cond = (res["r_start"] == -1)
    res.loc[cond, "r_start"] = res.loc[cond, cols].max(axis=1) + mock_insert_length
    
    cols = ["r_start", "r_end"]
    cond = (res["f_start"] == -1)
    res.loc[cond, "f_start"] = res.loc[cond, cols].min(axis=1) - mock_insert_length
    
    res["insert_length"] = res["r_start"] - res["f_start"]
    res["direction"] = res["insert_length"].apply(lambda x : "+" if x>0 else "-")
    res["Genes"] = ""
    
    for i, row in res.iterrows():
    
        if row["r_alignment_subject"] == row["f_alignment_subject"]:
            if (row["f_start"] > 0) & (row["r_start"] > 0):
                start = row["f_start"] if row["direction"] == "+" else row["r_start"]
                end = row["r_start"] if row["direction"] == "+" else row["f_start"]
                alignment_subject = row["f_alignment_subject"]
                
                gene_text = get_gene_text(df_annot, start, end, alignment_subject)
                res.loc[i, "Genes"] = gene_text
            else:
                res.loc[i, "Genes"] = "?"
                res.loc[i, "direction"] = "?"
                res.loc[i, "insert_length"] = "?"
            
            if (row["f_end"] < 0) | (row["r_end"] < 0):
                res.loc[i, "insert_length"] = "?"
                
        else:         
            # find the genes between start and end for forward
            alignment_subject = row['f_alignment_subject']
            start = min(row['f_start'], row['f_end'])
            end = max(row['f_start'], row['f_end'])
            gene_text_f = get_gene_text(df_annot, start, end, alignment_subject)
            
            # find the genes between start and end for reverse
            alignment_subject = row['r_alignment_subject']
            start = min(row['r_start'], row['r_end'])
            end = max(row['r_start'], row['r_end'])
            gene_text_r = get_gene_text(df_annot, start, end, alignment_subject)
            
            res.loc[i, 'Genes'] = gene_text_f + ", " + gene_text_r
            res.loc[i, 'direction'] = '?'
            res.loc[i, 'insert_length'] = '?'
            
    res.to_csv(outputfile, sep="\t", index=False)

    cols = ["Gene", "Name", "f_start", "f_end", "r_start", "r_end"]
    outputfile2 = outputfile.replace(".tsv", "_geneout.tsv")

    res["Gene"] = res["Genes"].str.split("\), ")
    geneout = res.explode("Gene")[cols].rename(columns={"Name": "Mapped_seguence"})
    geneout.to_csv(outputfile2, sep="\t", index=False)


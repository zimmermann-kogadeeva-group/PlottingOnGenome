#!/usr/bin/env python

from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from math import nan, isnan


def main(genome_file, annot_file, locus_file, output_prefix):
    # combined contigs downloaded from NCBI; only a single genome

    contigColors = {"blue": "#8DDEF7", "green": "#cffccc"}

    # Load geneMappingAlignment & geneLocations files
    geneMappingAlignment = pd.read_csv(annot_file, sep="\t")
    geneLocations = pd.read_csv(locus_file, sep=",")

    # replace "?" in "insert_length" column with "nan" so that it is also a floating point not a str. Important for later analysis.
    geneMappingAlignment["insert_length"] = (
        geneMappingAlignment["insert_length"]
        .str.replace("?", "nan", regex=False)
        .astype(float)
    )
    # delete all inserts that are larger than 1e6: nonsense inserts.
    geneMappingAlignment = geneMappingAlignment[
        np.absolute(geneMappingAlignment["insert_length"]) < 1e6
    ].reset_index()

    # create one list with all genome contigs (record id + record sequence) that were downloaded from NCBI
    genomeContigNames = []
    genomeContigLength = []
    for record in SeqIO.parse(genome_file, "fasta"):
        genomeContigNames.append(record.id)
        genomeContigLength.append(len(record.seq))

    # list where sequences were aligned to for fwd and rev primer read-out
    mappedHitContigs = list(geneMappingAlignment["f_alignment_subject"])
    mappedHitContigs.extend(list(geneMappingAlignment["r_alignment_subject"]))
    # leave only unique mapped contigs
    mappedHitContigs_dict = dict(
        (k, j) for j, k in enumerate(mappedHitContigs)
    )  # create dictionary
    idx_all = [j for j in list(mappedHitContigs_dict.values())]
    # now idx_all contain indeces of unique gene names
    mappedHitContigs = [mappedHitContigs[j] for j in idx_all]
    mappedHitContigs.sort()

    # sort contig names
    ind_dict = dict((k, i) for i, k in enumerate(genomeContigNames))
    inter = genomeContigNames[:]
    inter.sort()
    genomeContigNames = [genomeContigNames[ind_dict[x]] for x in inter]
    genomeContigLength = [genomeContigLength[ind_dict[x]] for x in inter]
    del inter, ind_dict

    # find genome contigs that contain inserts
    ind_dict = dict((k, i) for i, k in enumerate(genomeContigNames))
    # only keep contig names that are present in both lists (.intersection())
    inter = set(genomeContigNames).intersection(mappedHitContigs)
    inter = list(inter)
    inter.sort()
    indices = [ind_dict[x] for x in inter]

    # map sequences to contigs
    fig, ax = plt.subplots(figsize=(6, 30))

    # plot contigs and label only those that got mapping
    contigMaps = []
    contigStart = 0
    contigEnd = 0
    contigLength = 0
    genomeContigStart = []
    genomeContigEnd = []
    for i in range(len(genomeContigLength)):
        contigStart = contigEnd + 1
        contigLength = genomeContigLength[i]
        contigEnd = contigStart + contigLength
        if genomeContigNames[i] in inter:  # label only contigs that got mapping
            contigMaps.append(
                GraphicFeature(
                    start=contigStart,
                    end=contigEnd,
                    strand=+1,
                    color=contigColors["blue"],
                    label=genomeContigNames[i],
                )
            )
        else:
            contigMaps.append(
                GraphicFeature(
                    start=contigStart,
                    end=contigEnd,
                    strand=+1,
                    color=contigColors["blue"],
                )
            )

        genomeContigStart.append(
            contigStart
        )  # save contig start and end in the "whole genome"
        genomeContigEnd.append(contigEnd)

    for i in range(len(geneMappingAlignment)):
        f_contig = geneMappingAlignment["f_alignment_subject"][i]
        r_contig = geneMappingAlignment["r_alignment_subject"][i]
        # start is the start of the contig + coordinate in the contig
        mappingStart = (
            genomeContigStart[genomeContigNames.index(f_contig)]
            + geneMappingAlignment["f_start"][i]
        )
        mappingEnd = (
            genomeContigStart[genomeContigNames.index(f_contig)]
            + geneMappingAlignment["r_start"][i]
        )
        # only map if length is <10000
        if abs(mappingStart - mappingEnd) < 10000:
            contigMaps.append(
                GraphicFeature(
                    start=mappingStart,
                    end=mappingEnd,
                    strand=+1,
                    color=contigColors["green"],
                    label=geneMappingAlignment["Name"][i],
                )
            )
    # plot contigs
    record = CircularGraphicRecord(sequence_length=contigEnd, features=contigMaps)
    record.plot(ax=ax)

    # record.plot(fontsize = 8)
    fig.savefig(output_prefix + "_mapped_to_genome.pdf", transparent=True)
    fig.savefig(output_prefix + "_mapped_to_genome.png", transparent=True)

    #### PLOTTING ####
    # plot sequences onto the referene genome and display genes from they contain

    # extract gene names from gene locations file. Exclude "na" and convert it into a list
    geneLocationsNames = geneLocations["locus_tag"].dropna()

    # create new columns in the geneMappingAlignment array
    geneMappingAlignment["seq_start"] = [
        x["f_start"] if x["direction"] == "+" else x["r_start"]
        for i, x in geneMappingAlignment.iterrows()
    ]  # true start of the sequence
    geneMappingAlignment["seq_end"] = [
        x["r_start"] if x["direction"] == "+" else x["f_start"]
        for i, x in geneMappingAlignment.iterrows()
    ]  # true end of the sequence
    geneMappingAlignment["new_start"] = (
        geneMappingAlignment["seq_start"] - 5500
    )  # modified start of the sequence so that we can plot the area on the genome around the sequence
    geneMappingAlignment["new_end"] = (
        geneMappingAlignment["seq_end"] + 5500
    )  # modified end of the sequence so that we can plot the area on the genome around the sequence
    geneMappingAlignment["seg_length"] = (
        geneMappingAlignment["new_end"] - geneMappingAlignment["new_start"]
    )

    # extract on which contig your sequenced Fwd & Rev sequences are mapping to. Create separate lists for Fwd and Rev.
    f_contigNames = list(geneMappingAlignment["f_alignment_subject"])
    r_contigNames = list(geneMappingAlignment["r_alignment_subject"])

    # create pdf file and figure
    pdf = PdfPages(output_prefix + "_seqs_mapped_to_genome.pdf")

    # List all contig names that have sequences mapped to them
    # loop thorugh all mappedHitContigs and plot sequences located there onto the genome.
    for val in mappedHitContigs:
        if val in f_contigNames:
            curcontigNames = f_contigNames[:]
        else:
            curcontigNames = r_contigNames[:]

        geneMappingSubset = geneMappingAlignment[np.array(curcontigNames) == val]
        # split gene names to delete strand and percentage information. Only keep gene name.
        genes_split = (
            geneMappingSubset["Genes"]
            .str.split("\), ")
            .explode()
            .apply(lambda x: x.split(" ")[0] if type(x) == str else "")
        )

        geneLocationsSubset = geneLocations[
            geneLocations["locus_tag"].isin(genes_split)
        ]

        # sort rows according to "new_start" value so that sequences located at the beginning of the genome are displayed before those further downstream.
        geneMappingSubset.sort_values(by="new_start", inplace=True)

        # only display genes that are within the segment
        for i, row in geneMappingSubset.iterrows():
            # either the start or the end of a gene should be inside the sequence
            cond1 = (geneLocationsSubset["start"] >= row["seq_start"]) | (
                geneLocationsSubset["end"] >= row["seq_start"]
            )
            cond2 = (geneLocationsSubset["end"] <= row["seq_end"]) | (
                geneLocationsSubset["start"] <= row["seq_end"]
            )
            new_subset = geneLocationsSubset[cond1 & cond2]

            contigMappingFragments = []
            for j, row2 in new_subset.iterrows():
                if row2["start"] < row2["end"]:
                    contigMappingFragments.append(
                        GraphicFeature(
                            start=row2["start"],
                            end=row2["end"],
                            strand=+1,
                            color=contigColors["green"],
                            label=row2["locus_tag"],
                        )
                    )
                else:
                    contigMappingFragments.append(
                        GraphicFeature(
                            start=row2["end"],
                            end=row2["start"],
                            strand=-1,
                            color=contigColors["green"],
                            label=row2["locus_tag"],
                        )
                    )

            fig, axs = plt.subplots(2, 1, figsize=(8, 6))

            record = GraphicRecord(
                sequence_length=row["seg_length"],
                first_index=row["new_start"],
                features=contigMappingFragments,
            )
            record.plot(ax=axs[1], annotate_inline=True)

            # only sequences present in the seg_length should be displayed. No super large inserts spanning more than one segment should be displayed ("< 2*seg_length")
            cond1 = (geneMappingSubset["f_start"] >= row["new_start"]) | (
                geneMappingSubset["r_start"] >= row["new_start"]
            )
            cond2 = (geneMappingSubset["r_start"] <= row["new_end"]) | (
                geneMappingSubset["f_start"] <= row["new_end"]
            )
            cond3 = geneMappingSubset["seg_length"] < 2 * row["seg_length"]
            new_subset = geneMappingSubset[cond1 & cond2 & cond3]

            contigMappingFragments = []
            for j, row2 in new_subset.iterrows():
                if row2["f_start"] < row2["r_start"]:
                    contigMappingFragments.append(
                        GraphicFeature(
                            start=row2["f_start"],
                            end=row2["r_start"],
                            strand=+1,
                            color=contigColors["green"],
                            label=row2["Name"],
                        )
                    )
                else:
                    contigMappingFragments.append(
                        GraphicFeature(
                            start=row2["r_start"],
                            end=row2["f_start"],
                            strand=-1,
                            color=contigColors["green"],
                            label=row2["Name"],
                        )
                    )

            # plot contigs
            record = GraphicRecord(
                sequence_length=row["seg_length"],
                first_index=row["new_start"],
                features=contigMappingFragments,
            )
            record.plot(ax=axs[0], annotate_inline=True)

            axs[0].set_title(row["Name"])
            # set axis ticks to a max number of 8 ticks so that no overlaps in labels occur
            axs[0].xaxis.set_major_locator(plt.MaxNLocator(8))
            axs[1].xaxis.set_major_locator(plt.MaxNLocator(8))

            # save images as separate .png
            plt.savefig(
                output_prefix + row["Name"] + "_mapped_to_genome" + ".png",
                transparent=True,
            )
            # save appended images as .pdf
            pdf.savefig(fig)
    pdf.close()

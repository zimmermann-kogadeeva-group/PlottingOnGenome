from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from dna_features_viewer import (
    BiopythonTranslator,
    CircularGraphicRecord,
    GraphicFeature,
    GraphicRecord,
)


def set_feature(feature, **kwargs):
    """Helper function to set an attribute in Biopython feature"""
    new_feature = deepcopy(feature)
    for k, v in kwargs.items():
        if hasattr(new_feature, k):
            setattr(new_feature, k, v)
    return new_feature


def plot_insert(
    insert, genome, buffer=4000, figsize=None, axs=None, col1="#8DDEF7", col2="#CFFCCC"
):
    # Default values for figure size and create the figure
    if axs is None:
        figsize = figsize or (10, 8)
        fig, axs = plt.subplots(2, 1, figsize=figsize)
    assert len(axs) == 2

    # Create a new graphic object for query sequence
    features = [
        GraphicFeature(
            start=insert.start,
            end=insert.end,
            strand=insert.strand,
            color=col2,
            label=insert.seq_id,
        )
    ]

    # Plot the query sequence on the upper axes
    record_seq = GraphicRecord(
        first_index=insert.start - buffer,
        sequence_length=insert.end - insert.start + 2 * buffer,
        features=features,
    )
    _ = record_seq.plot(ax=axs[0])

    # Create graphic objects for all the genes and CDSes using
    # dna_features_viewer.BiopythonTranslator()
    conv = BiopythonTranslator()
    conv.default_feature_color = col1
    features = [
        conv.translate_feature(shift_feature(x, insert.start - buffer))
        for x in genome[insert.hit_id][
            insert.start - buffer : insert.end + buffer
        ].features
    ]

    # Plot the genes and CDSes in the region of the mapped sequence
    record_hits = GraphicRecord(
        first_index=insert.start - buffer,
        sequence_length=insert.end - insert.start + 2 * buffer,
        features=features,
    )
    _ = record_hits.plot(ax=axs[1])

    return axs


def plot_on_genome(
    inserts,
    genome,
    labels=True,
    figsize=None,
    col1="#8DDEF7",
    col2="#CFFCCC",
    ax=None,
):
    if ax is None:
        figsize = figsize or (10, 20)
        fig, ax = plt.subplots(figsize=figsize)

    # Get just the sequences for each NCBI record and order them by size in
    # descending order. 'x.features[0]' to get the whole sequence for a
    # given NCBI record. Other features are specific genes, cfs, etc.
    db_seqs = sorted(
        [x for x in genome.values() if x.seq.defined],
        key=lambda x: len(x),
        reverse=True,
    )

    # Get the shifts needed to plot all the NCBI records in a continuous line
    shifts = np.cumsum([0] + [len(x) for x in db_seqs])

    # Get IDs of NCBI records that were mapped to. Used to check where to
    # add labels if option is set.
    mapped_ids = {x.hit_id for x in inserts}

    # Make plots of NCBI records and label only the ones that were mapped
    # to. Using BiopythonTranslator() didn't allow for control of labels,
    # hence we are just using GraphicFeature class
    features = [
        (
            GraphicFeature(
                start=shifts[i],
                end=shifts[i + 1],
                label=x.id if (labels and x.id in mapped_ids) else None,
                color=col1,
            )
        )
        for i, x in enumerate(db_seqs)
    ]

    # Get IDs of NCBI records in the order as in the figure. Used to make
    # sure locations are shifted correctly.
    ids = [x.id for x in db_seqs]

    # Add plots of the query sequences plotted on top of the plots of NCBI records
    hits = [
        GraphicFeature(
            start=insert.start + shifts[ids.index(insert.hit_id)] + 1,
            end=insert.end + shifts[ids.index(insert.hit_id)],
            strand=insert.strand,
            color=col2,
            label=None,
        )
        for insert in inserts
    ]

    rec = CircularGraphicRecord(sequence_length=shifts[-1], features=features + hits)

    _ = rec.plot(ax, annotate_inline=False)

    return ax


def plot_insert_dists(inserts, axs=None):
    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 3, figsize=(12, 5))
    assert len(axs) == 3

    insert_lengths = [len(x) for x in inserts]
    props = [x for x in inserts for x in x.coverage]

    if insert_lengths:
        sns.histplot(x=insert_lengths, ax=axs[0])
        axs[0].set(title="Insert lengths")

        sns.stripplot(insert_lengths, ax=axs[1], color="black")
        sns.violinplot(insert_lengths, ax=axs[1], width=0.5, saturation=0.4, inner=None)
        sns.boxplot(insert_lengths, ax=axs[1], width=0.25, boxprops={"zorder": 2})
        axs[1].set(xticklabels=[], title="Insert lengths")

    if props:
        sns.stripplot(props, ax=axs[2], color="black")
        sns.violinplot(props, ax=axs[2], width=0.5, saturation=0.4, inner=None)
        sns.boxplot(props, ax=axs[2], width=0.2, boxprops={"zorder": 2})
        axs[2].set(xticklabels=[], title="Proportions")

    return axs

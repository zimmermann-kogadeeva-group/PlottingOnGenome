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


def _get_graphic_records_insert(insert, buffer, col1, col2):
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

    # Create graphic objects for all the genes and CDSes using
    # dna_features_viewer.BiopythonTranslator()
    conv = BiopythonTranslator()
    conv.default_feature_color = col1
    features = [conv.translate_feature(x) for x in insert.genes]

    # Plot the genes and CDSes in the region of the mapped sequence
    record_hits = GraphicRecord(
        first_index=insert.start - buffer,
        sequence_length=insert.end - insert.start + 2 * buffer,
        features=features,
    )

    return record_seq, record_hits


def plot_insert(
    insert,
    buffer=4000,
    col1="#8DDEF7",
    col2="#CFFCCC",
    figsize=None,
    axs=None,
    backend="matplotlib",
):

    seqs, hits = _get_graphic_records_insert(insert, buffer, col1, col2)

    figsize = figsize or (10, 8)

    if backend == "matplotlib":
        # Default values for figure size and create the figure
        if axs is None:
            fig, axs = plt.subplots(2, 1, figsize=figsize)
        assert len(axs) == 2

        # Create a new graphic object for query sequence
        _ = seqs.plot(ax=axs[0])
        _ = hits.plot(ax=axs[1])

        return axs


def _get_contig_label(contig, mapped_ids, show_labels=True):
    if show_labels:
        if mapped_ids and contig.id in mapped_ids:
            return contig.id
        else:
            return None
    else:
        return None


def _get_graphic_records_genome(genome, inserts, show_labels, col1, col2):

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
                label=_get_contig_label(x, mapped_ids, show_labels),
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

    return rec


def plot_genome(
    genome,
    inserts=None,
    show_labels=True,
    col1="#8DDEF7",
    col2="#CFFCCC",
    figsize=None,
    ax=None,
    backend="matplotlib",
):
    rec = _get_graphic_records_genome(genome, inserts, show_labels, col1, col2)

    figsize = figsize or (10, 20)

    if backend == "matplotlib":
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        inserts = inserts or []

        _ = rec.plot(ax, annotate_inline=False)

    return ax


def plot_histogram(inserts, axs=None):
    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    assert len(axs) == 2

    insert_lengths = [len(x) for x in inserts]
    coverage = [x for x in inserts for x in x.coverage]

    if insert_lengths:
        sns.histplot(x=insert_lengths, ax=axs[0])
        axs[0].set(title="Insert lengths")
    if coverage:
        sns.histplot(x=coverage, ax=axs[1])
        axs[1].set(title="Coverage")


def plot_dists(inserts, color="steelblue", saturation=0.4, width=1.0, axs=None):
    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    assert len(axs) == 2

    insert_lengths = [len(x) for x in inserts]
    coverage = [x for x in inserts for x in x.coverage]

    if insert_lengths:
        sns.stripplot(x=1, y=insert_lengths, ax=axs[0], color="black")
        sns.violinplot(
            x=2,
            y=insert_lengths,
            ax=axs[0],
            color=color,
            width=width,
            saturation=saturation,
            inner=None,
        )
        sns.boxplot(
            x=3,
            y=insert_lengths,
            ax=axs[0],
            color=color,
            width=0.5 * width,
            boxprops={"zorder": 2},
        )
        axs[0].set(xticklabels=[], title="Insert lengths")

    if coverage:
        sns.stripplot(x=1, y=coverage, ax=axs[1], color="black")
        sns.violinplot(
            x=2,
            y=coverage,
            ax=axs[1],
            color=color,
            width=width,
            saturation=saturation,
            inner=None,
        )
        sns.boxplot(
            x=3,
            y=coverage,
            ax=axs[1],
            color=color,
            width=0.5 * width,
            boxprops={"zorder": 2},
        )
        axs[1].set(xticklabels=[], title="Coverage")

    return axs

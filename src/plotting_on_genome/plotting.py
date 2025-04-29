from copy import deepcopy
from itertools import chain, cycle

import matplotlib.pyplot as plt
import seaborn as sns
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


def set_feature(feature, **kwargs):
    """Helper function to set an attribute in Biopython feature"""
    new_feature = deepcopy(feature)
    for k, v in kwargs.items():
        if hasattr(new_feature, k):
            setattr(new_feature, k, v)
    return new_feature


def fig_axvline(axes, value, ls="--", color="gray", **kwargs):
    fig = axes[0].get_figure()
    transFigure = fig.transFigure.inverted()

    coord1 = transFigure.transform(axes[0].transData.transform([value, 0]))
    coord2 = transFigure.transform(axes[1].transData.transform([value, 0]))

    line = plt.Line2D(
        (coord1[0], coord2[0]),
        (coord1[1], coord2[1]),
        ls=ls,
        color=color,
        transform=fig.transFigure,
        zorder=-100,
    )
    fig.lines.append(line)


def plot_histogram(inserts, axs=None, **kwargs):
    if "figsize" not in kwargs:
        kwargs["figsize"] = (12, 5)

    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 2, **kwargs)
    assert len(axs) == 2

    insert_lengths = [len(x) for x in inserts]
    coverage = [x.coverage for x in inserts]

    if insert_lengths:
        sns.histplot(x=insert_lengths, ax=axs[0])
        axs[0].set(title="Insert lengths")
    if coverage:
        sns.histplot(x=coverage, ax=axs[1])
        axs[1].set(title="Coverage")


def plot_dists(
    inserts, color="steelblue", saturation=0.4, width=1.0, axs=None, **kwargs
):
    if "figsize" not in kwargs:
        kwargs["figsize"] = (12, 5)

    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 2, **kwargs)
    assert len(axs) == 2

    insert_lengths = [len(x) for x in inserts]
    coverage = [x.coverage for x in inserts]

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


def plot_multiple_genomes(
    *inserts_dict,
    seq_id=None,
    insert_type="both",
    filter_threshold=None,
    palette="tab10",
    ax=None,
    **kwargs,
):
    if "figsize" not in kwargs:
        kwargs["figsize"] = (10, 8)

    if ax is None:
        fig, ax = plt.subplots(**kwargs)

    palette = color_sequences[palette]
    message = "Number of colors in palette is not sufficient for number of genomes"
    assert len(inserts_dict) <= len(palette), message

    res = [
        ins_dict.get_graphic_features(
            seq_id, insert_type, filter_threshold, col1=col, col2=col, linecolor=col
        )
        for ins_dict, col in zip(inserts_dict, palette)
    ]

    features = list(chain.from_iterable([x[1] for x in res]))
    rec = CircularGraphicRecord(max([x[0] for x in res]), features)

    _ = rec.plot(ax, annotate_inline=False)

    return ax

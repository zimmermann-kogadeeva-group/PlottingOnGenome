from copy import deepcopy

import matplotlib.pyplot as plt
import seaborn as sns


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


def plot_histogram(inserts, axs=None):
    # Default values for figure size and create the figure
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    assert len(axs) == 2

    insert_lengths = [len(x) for x in inserts]
    coverage = [x.coverage for x in inserts]

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

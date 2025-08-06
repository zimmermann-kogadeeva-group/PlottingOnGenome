from collections import defaultdict
from itertools import chain, cycle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


def clusters_to_seq_labels(clusters):
    # Get sequences / insert labels based on clusters
    seq_labels = defaultdict(dict)
    if clusters is not None:
        for (genome, clust_idx), seq_ids in clusters.items():
            for seq_id in seq_ids:
                seq_labels[genome][seq_id] = f"Cluster {clust_idx}"
    return seq_labels


def merge_sets(set1, set2):
    set1_cond = set1 is None or not len(set1)
    set2_cond = set2 is None or not len(set2)
    if set1_cond and not set2_cond:
        return set2
    elif set2_cond and not set1_cond:
        return set1
    elif set1_cond and set2_cond:
        return None
    else:
        return list(set([*set1, *set2]))


def merge_selections(selection, seq_labels):
    keys = set([*selection.keys(), *seq_labels.keys()])
    return {k: merge_sets(selection[k], list(seq_labels[k].keys())) for k in keys}


class Comparison(dict):

    def __init__(self, *args, **kwargs):
        super(Comparison, self).__init__(*args, **kwargs)

    def get_inserts(
        self, selection=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if selection is None:
            selection = list(self.keys())

        if isinstance(selection, (list, tuple)):
            selection = {x: None for x in selection}

        return [
            self.__getitem__(sel).get(
                insert_ids,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
            )
            for sel, insert_ids in selection.items()
        ]

    def get_inserts_df(
        self, selection=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if selection is None:
            selection = self.keys()

        inserts_dfs = [
            self.__getitem__(sel)
            .to_dataframe(insert_type=insert_type, filter_threshold=filter_threshold)
            .assign(genome=sel)
            for sel in selection
        ]

        if len(inserts_dfs):
            return pd.concat(inserts_dfs, ignore_index=True)

    def get_insert_presence_df(
        self, selection=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if selection is None:
            selection = self.keys()

        dfs = [
            pd.DataFrame(
                self.__getitem__(sel).get_seq_ids(insert_type, filter_threshold),
                columns=["insert_ids"],
            ).assign(genome=sel)
            for sel in selection
        ]

        df_insert_presence = (
            pd.concat(dfs, ignore_index=True)
            .groupby(["insert_ids", "genome"], as_index=False)
            .agg(num_inserts=pd.NamedAgg("insert_ids", "count"))
            .pivot(index="insert_ids", columns="genome", values="num_inserts")
        )

        return df_insert_presence

    def get_clusters(
        self, selection=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if selection is None:
            selection = self.keys()

        return {
            (sel, idx): cluster
            for sel in selection
            for idx, cluster in enumerate(
                self.__getitem__(sel).get_clusters(insert_type, filter_threshold)
            )
        }

    def get_genes_df(
        self,
        keys=None,
        seq_id=None,
        insert_type="both",
        filter_threshold=None,
        buffer=4000,
        **kwargs,
    ):
        if keys is None:
            keys = self.keys()

        genes_dfs = [
            self.__getitem__(key)
            .genes_to_dataframe(
                seq_id,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                buffer=buffer,
            )
            .assign(genome=key)
            for key in keys
        ]

        if len(genes_dfs):
            return pd.concat(genes_dfs, ignore_index=True)

    def plot(
        self,
        selection=None,
        clusters=None,
        insert_type="both",
        filter_threshold=None,
        show_labels=True,
        show_titles=False,
        palette="tab10",
        facet_wrap=None,
        genomes_order=None,
        fig=None,
        **kwargs,
    ):
        if selection is None:
            selection = list(self.keys())

        if isinstance(selection, (list, tuple)):
            selection = {x: None for x in selection}

        if genomes_order is None:
            genomes_order = list(selection.keys())

        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        num_genomes = len(selection)
        nrows, ncols = 1, 1
        if facet_wrap is not None:
            ncols = min(facet_wrap, num_genomes)
            nrows = num_genomes // ncols

        if fig is None:
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        else:
            axs = fig.get_axes()

        # Check that axs is iterable even if a single axes object
        if isinstance(axs, plt.Axes):
            axs = np.array([axs])

        palette = cycle(color_sequences[palette])
        seq_labels = clusters_to_seq_labels(clusters)

        # update selection based on clusters
        selection = merge_selections(selection, seq_labels)

        res = [
            self.__getitem__(key).get_graphic_features(
                selection[key],
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                show_labels=show_labels,
                seq_labels=seq_labels[key],
                col1=col,
                col2=col,
                linecolor=col,
            )
            for key, col in zip(genomes_order, palette)
        ]

        genomes_labels = genomes_order
        if facet_wrap is None:
            features = list(chain.from_iterable([x[1] for x in res]))
            rec = CircularGraphicRecord(max([x[0] for x in res]), features)
            _ = rec.plot(axs[0], annotate_inline=False)
            genomes_labels = [", ".join(genomes_order)]

        else:
            for (seq_len, features), ax in zip(res, axs.flatten()):
                rec = CircularGraphicRecord(seq_len, features)
                _ = rec.plot(ax, annotate_inline=False)

            # Remove redundant axis
            for ax in axs[num_genomes:]:
                ax.axis("off")

        # Set titles
        if show_titles:
            for title, ax in zip(genomes_labels, axs):
                ax.set_title(title)

        return fig

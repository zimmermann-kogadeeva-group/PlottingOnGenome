from itertools import chain, cycle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


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

    # TODO: This needs to be rethought
    def get_inserts(
        self,
        genomes=None,
        seq_ids=None,
        insert_ids=None,
        insert_type="both",
        filter_threshold=None,
        **kwargs,
    ):
        if genomes is None:
            genomes = self.keys()

        return [
            self.__getitem__(g).get(
                insert_ids=insert_ids,
                seq_ids=seq_ids,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
            )
            for g in genomes
        ]

    def get_inserts_df(
        self, genomes=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if genomes is None:
            genomes = self.keys()

        inserts_dfs = [
            self.__getitem__(genome)
            .to_dataframe(insert_type=insert_type, filter_threshold=filter_threshold)
            .assign(genome=genome)
            for genome in genomes
        ]

        if len(inserts_dfs):
            return pd.concat(inserts_dfs, ignore_index=True)

    def get_insert_presence_df(
        self, genomes=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if genomes is None:
            genomes = self.keys()

        dfs = [
            pd.DataFrame(
                self.__getitem__(g).get_seq_ids(insert_type, filter_threshold),
                columns=["insert_ids"],
            ).assign(genome=g)
            for g in genomes
        ]

        df_insert_presence = (
            pd.concat(dfs, ignore_index=True)
            .groupby(["insert_ids", "genome"], as_index=False)
            .agg(num_inserts=pd.NamedAgg("insert_ids", "count"))
            .pivot(index="insert_ids", columns="genome", values="num_inserts")
        )

        return df_insert_presence

    # TODO: this needs to be rethought as well (structure of the output)
    def get_clusters(
        self, genomes=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if genomes is None:
            genomes = self.keys()

        return {
            g: self.__getitem__(g).get_clusters(insert_type, filter_threshold)
            for g in genomes
        }

    # TODO: this needs to be rethought as well (structure of the output)
    def get_genes_df(
        self,
        seq_ids=None,
        insert_ids=None,
        insert_type="both",
        filter_threshold=None,
        buffer=4000,
        **kwargs,
    ):
        if seq_ids is None:
            seq_ids = {x: None for x in self.keys()}
        elif isinstance(seq_ids, (list, tuple)):
            seq_ids = {x: seq_ids for x in self.keys()}

        if insert_ids is None:
            insert_ids = {x: None for x in self.keys()}
        elif isinstance(insert_ids, (list, tuple)):
            insert_ids = {x: insert_ids for x in self.keys()}

        genomes = set(seq_ids.keys()) | set(insert_ids.keys())

        genes_dfs = [
            self.__getitem__(g)
            .genes_to_dataframe(
                seq_ids=seq_ids.get(g),
                insert_ids=insert_ids.get(g),
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                buffer=buffer,
            )
            .assign(genome=g)
            for g in genomes
        ]

        if len(genes_dfs):
            drop_by = ["start", "end", "type", "insert_idx", "seq_id", "genome"]
            return pd.concat(genes_dfs, ignore_index=True).drop_duplicates(
                subset=drop_by
            )

    def get_graphic_features(
        self,
        seq_ids=None,
        insert_ids=None,
        genomes_order=None,
        seq_labels=None,
        contig_labels=True,
        insert_type="both",
        filter_threshold=None,
        palette="tab10",
    ):
        if seq_ids is None:
            seq_ids = {x: None for x in self.keys()}
        elif isinstance(seq_ids, (list, tuple)):
            seq_ids = {x: seq_ids for x in self.keys()}

        if insert_ids is None:
            insert_ids = {x: None for x in self.keys()}
        elif isinstance(insert_ids, (list, tuple)):
            insert_ids = {x: insert_ids for x in self.keys()}

        if genomes_order is None:
            genomes_order = set(seq_ids.keys()) | set(insert_ids.keys())

        if seq_labels is None:
            seq_labels = {x: None for x in self.keys()}

        palette = cycle(color_sequences[palette])

        res = {
            g: self.__getitem__(g).get_graphic_features(
                seq_ids=seq_ids.get(g),
                insert_ids=insert_ids.get(g),
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                contig_labels=contig_labels,
                seq_labels=seq_labels.get(g),
                col1=col,
                col2=col,
                linecolor=col,
            )
            for g, col in zip(genomes_order, palette)
        }

        return res

    # TODO: this needs to be rethought as well (structure of the output)
    def plot(
        self,
        seq_ids=None,
        insert_ids=None,
        genomes_order=None,
        seq_labels=None,
        contig_labels=True,
        insert_type="both",
        filter_threshold=None,
        show_titles=False,
        palette="tab10",
        facet_wrap=None,
        fig=None,
        **kwargs,
    ):
        # Get the graphic feature objects
        res = self.get_graphic_features(
            seq_ids,
            insert_ids,
            genomes_order,
            seq_labels,
            contig_labels,
            insert_type,
            filter_threshold,
            palette,
        )

        # Set all the plotting options based on input
        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        genome_labels = res.keys()
        num_genomes = len(genome_labels)

        # Set the number of axes in the figure based on facet_wrap option
        nrows, ncols = 1, 1
        if facet_wrap is not None:
            ncols = min(facet_wrap, num_genomes)
            nrows = num_genomes // ncols

        # Create figure and axes if not provided
        if fig is None:
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
        else:
            axs = fig.get_axes()

        # Check that axs is iterable even if a single axes object
        if isinstance(axs, plt.Axes):
            axs = np.array([axs])

        # Make and label all the plots
        if facet_wrap is None:
            features = list(chain.from_iterable([x[1] for x in res.values()]))
            rec = CircularGraphicRecord(max([x[0] for x in res.values()]), features)
            _ = rec.plot(axs[0], annotate_inline=False)
            axs[0].set_title(", ".join(list(genome_labels)))

        else:
            for (label, (seq_len, features)), ax in zip(res.items(), axs.flatten()):
                rec = CircularGraphicRecord(seq_len, features)
                _ = rec.plot(ax, annotate_inline=False)
                ax.set_title(label)

            # Remove redundant axis
            for ax in axs[num_genomes:]:
                ax.axis("off")

        return fig

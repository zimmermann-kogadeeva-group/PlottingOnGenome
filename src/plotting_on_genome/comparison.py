from collections import defaultdict
from itertools import chain, cycle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


class Clusters(dict):
    def __init__(self, clusters, genomes=None):
        super(Clusters, self).__init__(clusters)

        self.genomes = genomes or {g for g, clust_idx in self.keys()}
        self.cluster_labels = list(self.keys())

    def __contains__(self, item):
        return item in self.cluster_labels

    @property
    def insert_ids(self):
        insert_ids = {
            genome: list(
                chain.from_iterable(
                    (
                        insert_ids
                        for (g, clust_idx), insert_ids in self.items()
                        if g == genome
                    )
                )
            )
            for genome in self.genomes
        }
        return insert_ids

    @property
    def insert_labels(self):
        seq_labels = defaultdict(dict)
        for g, clust_idx in self.cluster_labels:
            insert_ids = self.__getitem__((g, clust_idx))
            # Set it such that only one label is given per cluster
            for insert_id in insert_ids:
                seq_labels[g][insert_id] = None
            seq_labels[g][insert_ids[-1]] = f"Cluster {clust_idx}"
        return seq_labels

    @property
    def other_repr(self):
        other = defaultdict(list)
        for (g, clust_idx), insert_ids in self.items():
            other[g].append((clust_idx, insert_ids))
        return other

    def subset(self, selection):
        return Clusters(
            {
                (g, clust_idx): self.__getitem__((g, clust_idx))
                for g, clust_idx in selection
            }
        )


class Comparison(dict):

    def __init__(self, *args, **kwargs):
        super(Comparison, self).__init__(*args, **kwargs)

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

        all_seq_ids = sorted(
            set(chain.from_iterable(self.__getitem__(g).seq_ids for g in genomes))
        )

        df_insert_presence = (
            pd.concat(dfs, ignore_index=True)
            .groupby(["insert_ids", "genome"], as_index=False)
            .agg(num_inserts=pd.NamedAgg("insert_ids", "count"))
            .pivot(index="insert_ids", columns="genome", values="num_inserts")
            .reindex(columns=genomes, index=all_seq_ids)
        )

        return df_insert_presence

    def get_clusters(
        self,
        genomes=None,
        insert_type="both",
        filter_threshold=None,
        plain_dict=True,
        **kwargs,
    ):
        if genomes is None:
            genomes = self.keys()

        if plain_dict:
            return {
                g: self.__getitem__(g).get_clusters(insert_type, filter_threshold)
                for g in genomes
            }
        else:
            return Clusters(
                {
                    (g, clust_idx): cluster
                    for g in genomes
                    for clust_idx, cluster in enumerate(
                        self.__getitem__(g).get_clusters(insert_type, filter_threshold)
                    )
                },
                genomes=genomes,
            )

    def get_labels(
        self,
        seq_ids,
        insert_ids=None,
        genomes=None,
        insert_type="both",
        filter_threshold=None,
    ):
        if genomes is None:
            genomes = self.keys()

        if insert_ids is None:
            insert_ids = {g: None for g in genomes}

        if seq_ids is not None or insert_ids is not None:
            return {
                g: self.__getitem__(g).get_labels(
                    seq_ids, insert_ids.get(g), genomes, insert_type, filter_threshold
                )
                for g in genomes
            }

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
            df_genes = pd.concat(genes_dfs, ignore_index=True).drop_duplicates(
                subset=drop_by
            )
            df_genes.insert(0, "genome", df_genes.pop("genome"))
            return df_genes

    def get_graphic_features(
        self,
        seq_ids=None,
        insert_ids=None,
        genomes_order=None,
        seq_labels=None,
        show_contig_labels=True,
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
                contig_labels=show_contig_labels,
                seq_labels=seq_labels.get(g),
                col1=col,
                col2=col,
                linecolor=col,
            )
            for g, col in zip(genomes_order, palette)
        }

        return res

    def plot(
        self,
        seq_ids=None,
        insert_ids=None,
        genomes_order=None,
        seq_labels=None,
        insert_type="both",
        filter_threshold=None,
        show_contig_labels=True,
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
            show_contig_labels,
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

            if show_titles:
                axs[0].set_title(", ".join(list(genome_labels)))

        else:
            for (label, (seq_len, features)), ax in zip(res.items(), axs.flatten()):
                rec = CircularGraphicRecord(seq_len, features)
                _ = rec.plot(ax, annotate_inline=False)

                if show_titles:
                    ax.set_title(label)

            # Remove redundant axis
            for ax in axs[num_genomes:]:
                ax.axis("off")

        return fig

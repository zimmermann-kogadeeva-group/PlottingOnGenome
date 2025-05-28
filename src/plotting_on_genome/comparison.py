from itertools import chain, cycle

import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


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
        buffer=None,
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
        insert_type="both",
        filter_threshold=None,
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

        palette = cycle(color_sequences[palette])

        res = [
            self.__getitem__(key).get_graphic_features(
                selection[key],
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                col1=col,
                col2=col,
                linecolor=col,
            )
            for key, col in zip(genomes_order, palette)
        ]

        if facet_wrap is None:
            features = list(chain.from_iterable([x[1] for x in res]))
            rec = CircularGraphicRecord(max([x[0] for x in res]), features)
            _ = rec.plot(axs, annotate_inline=False)
        else:
            for (seq_len, features), ax in zip(res, axs.flatten()):
                rec = CircularGraphicRecord(seq_len, features)
                _ = rec.plot(ax, annotate_inline=False)

            # Remove redundant axis
            for ax in axs.flatten()[num_genomes:]:
                ax.axis("off")

        return fig

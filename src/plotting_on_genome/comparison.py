from itertools import chain

import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import CircularGraphicRecord
from matplotlib import color_sequences


class Comparison(dict):

    def __init__(self, *args, **kwargs):
        super(Comparison, self).__init__(*args, **kwargs)

    def get_inserts_df(
        self, keys=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if keys is None:
            keys = self.keys()

        inserts_dfs = [
            self.__getitem__(key)
            .to_dataframe(insert_type=insert_type, filter_threshold=filter_threshold)
            .assign(genome=key)
            for key in keys
        ]

        if len(inserts_dfs):
            return pd.concat(inserts_dfs, ignore_index=True)

    def get_insert_presence_df(
        self, keys=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        if keys is None:
            keys = self.keys()

        dfs = [
            pd.DataFrame(
                self.__getitem__(key).get_insert_ids(insert_type, filter_threshold),
                columns=["insert_ids"],
            ).assign(genome=key)
            for key in keys
        ]

        df_insert_presence = (
            pd.concat(dfs, ignore_index=True)
            .groupby(["insert_ids", "genome"], as_index=False)
            .agg(num_inserts=pd.NamedAgg("insert_ids", "count"))
            .pivot(index="insert_ids", columns="genome", values="num_inserts")
        )

        return df_insert_presence

    def get_genes_df(
        self,
        keys=None,
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
        keys=None,
        seq_id=None,
        insert_type="both",
        filter_threshold=None,
        palette="tab10",
        ax=None,
        **kwargs,
    ):
        if keys is None:
            keys = self.keys()

        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        if ax is None:
            fig, ax = plt.subplots(**kwargs)

        palette = color_sequences[palette]
        message = "Number of colors in palette is not sufficient for number of genomes"
        assert len(self) <= len(palette), message

        res = [
            self.__getitem__(key).get_graphic_features(
                seq_id, insert_type, filter_threshold, col1=col, col2=col, linecolor=col
            )
            for key, col in zip(keys, palette)
        ]

        features = list(chain.from_iterable([x[1] for x in res]))
        rec = CircularGraphicRecord(max([x[0] for x in res]), features)

        _ = rec.plot(ax, annotate_inline=False)

        return ax

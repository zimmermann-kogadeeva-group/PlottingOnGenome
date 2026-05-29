from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
)
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap

from .helper import (
    fig_axvline,
    get_gene_coverage,
    get_genes_graphic_record,
    shift_feature,
)


def set_feature(feature, **kwargs):
    """Helper function to set an attribute in Biopython feature"""
    new_feature = deepcopy(feature)
    for k, v in kwargs.items():
        if hasattr(new_feature, k):
            setattr(new_feature, k, v)
    return new_feature


def _get_endpoints(hsp1, hsp2=None, avg_insert_length=0):
    # TODO: check seq_len
    if hsp2 is not None:
        if hsp1.hit_strand == 1:
            start = hsp1.hit_start
            middle1 = hsp1.hit_end
            middle2 = hsp2.hit_start
            end = hsp2.hit_end
        else:
            start = hsp2.hit_start
            middle1 = hsp2.hit_end
            middle2 = hsp1.hit_start
            end = hsp1.hit_end
        return start, middle1, middle2, end
    else:
        start = hsp1.hit_start
        end = max(hsp1.hit_end, hsp1.hit_start + avg_insert_length)
        middle1 = middle2 = hsp1.hit_end
        return start, middle1, middle2, end


class Insert(object):

    def __init__(
        self,
        seq_id,
        idx,
        seq1,
        hsp1,
        seq2=None,
        hsp2=None,
        genome=None,
        avg_insert_len=4000,
        paired=True,
    ):
        # TODO: maybe switch to using a list of hsps instead of hsp1 and hsp2
        self.hsp1 = hsp1
        self.hsp2 = hsp2
        self.seq1 = seq1
        self.seq2 = seq2
        self.strand = self.hsp1.hit_strand
        self.hit_id = self.hsp1.hit_id
        self.query_id = self.hsp1.query_id
        self.seq_id = seq_id
        self.idx = idx
        self.paired = True if hsp2 is not None else False

        if paired:
            self.start, self.middle1, self.middle2, self.end = _get_endpoints(
                hsp1, hsp2, avg_insert_len
            )
        else:
            self.start, self.end = self.hsp1.hit_start, self.hsp1.hit_end
            self.middle1 = self.middle2 = None

        self._genome = genome

        self.cov1 = len(hsp1.query.seq) / len(seq1)
        self.cov2 = len(hsp2.query.seq) / len(seq2) if hsp2 is not None else None

        self.coverage = min(x for x in (self.cov1, self.cov2) if x is not None)

    def __repr__(self):
        return (
            "Insert("
            f"seq_id={self.seq_id}, "
            f"idx={self.idx}, "
            f"start={self.start}, "
            f"end={self.end}, "
            f"hit_id={self.hit_id}, "
            f"coverage={self.coverage:.2f}, "
            f"paired={self.paired})"
        )

    def __format__(self, format_spec):
        if format_spec == "short":
            return (
                f"seq_id = {self.seq_id}, "
                f"index = {self.idx}, "
                f"coverage = {self.coverage:.2f}, "
                f"paired = {self.paired}, "
                f"hit_id = {self.hit_id}"
            )
        return self.__repr__()

    def __len__(self):
        return self.end - self.start + 1

    def get_genes(self, buffer=4000):
        if isinstance(buffer, int):
            buffer = (buffer, buffer)
        assert isinstance(buffer, (tuple, list)) and len(buffer) == 2

        start = self.start - buffer[0]
        end = self.end + buffer[1]

        return [shift_feature(gene, start) for gene in self._genome[start:end].features]

    def _get_gene_coverage(self, gene):
        return get_gene_coverage(
            gene.location.start, gene.location.end, self.start, self.end
        )

    def to_dataframe(self, buffer=4000):
        if isinstance(buffer, int):
            buffer = (buffer, buffer)
        assert isinstance(buffer, (tuple, list)) and len(buffer) == 2

        return pd.DataFrame(
            [
                {
                    "start": gene.location.start,
                    "end": gene.location.end,
                    "strand": gene.location.strand,
                    "type": gene.type,
                    "coverage": self._get_gene_coverage(gene),
                    **gene.qualifiers,
                }
                for gene in self.get_genes(buffer)
            ]
        )

    def get_graphic_feature(self, col="#2e8b57"):
        return GraphicFeature(
            start=self.start,
            end=self.end,
            strand=self.strand,
            color=col,
            label=self.seq_id,
        )

    # TODO: col1, col2, cmap - simplify
    def get_genes_graphic_record(
        self,
        buffer=4000,
        col="#ebf3ed",
        feature_types=None,
        cmap=None,
    ):
        # Get genes close to the given insert
        genes = self.get_genes(buffer)

        # Get the GraphicRecord obj to be plotted
        return get_genes_graphic_record(
            genes,
            self.start,
            self.end,
            buffer=buffer,
            col=col,
            feature_types=feature_types,
            cmap=cmap,
        )

    def plot(
        self,
        buffer=4000,
        col1="#ebf3ed",
        col2="#2e8b57",
        feature_types=None,
        colorbar=False,
        show_mappings=False,
        axs=None,
        backend="matplotlib",
        **kwargs,
    ):

        if isinstance(buffer, int):
            buffer = (buffer, buffer)
        assert isinstance(buffer, (tuple, list)) and len(buffer) == 2

        cmap = None
        if colorbar:
            cmap = LinearSegmentedColormap.from_list("custom", [col1, col2])

        # Plot the query sequence on the upper axes
        rec_seq = GraphicRecord(
            first_index=self.start - buffer[0],
            sequence_length=self.end - self.start + sum(buffer),
            features=[self.get_graphic_feature(col2)],
        )
        rec_genes = self.get_genes_graphic_record(buffer, col1, feature_types, cmap)

        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        if backend == "matplotlib":
            # Default values for figure size and create the figure
            if axs is None:
                fig, axs = plt.subplots(2, 1, **kwargs)
            else:
                fig = axs[0].get_figure()
            assert len(axs) == 2

            # Create a new graphic object for query sequence
            _ = rec_seq.plot(ax=axs[0])
            _ = rec_genes.plot(ax=axs[1])

            if colorbar:
                ax_cb = fig.add_axes([0.9, 0.1, 0.02, 0.8])
                fig.subplots_adjust(left=0.1, right=0.85)
                ColorbarBase(
                    ax_cb, cmap=cmap, orientation="vertical", label="gene coverage"
                )

            fig_axvline(axs, self.start)
            fig_axvline(axs, self.end)
            if show_mappings:
                fig_axvline(axs, self.middle1)
                fig_axvline(axs, self.middle2)

            return axs
        else:
            raise NotImplementedError

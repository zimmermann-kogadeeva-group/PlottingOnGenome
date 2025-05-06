from copy import deepcopy

import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import (
    BiopythonTranslator,
    GraphicFeature,
    GraphicRecord,
)
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap

from .helper import shift_feature


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


class Insert(object):

    def _get_endpoints(self, hsp1, hsp2=None, avg_insert_length=0):
        # TODO: check seq_len
        if hsp2 is not None:
            if hsp1.hit_strand == 1:
                start = hsp1.hit_start
                end = hsp2.hit_end
            else:
                start = hsp2.hit_start
                end = hsp1.hit_end
            return start, end
        else:
            start = hsp1.hit_start
            end = hsp1.hit_start + avg_insert_length
            return start, end

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
        self.matched = True if hsp2 is not None else False

        if paired:
            self.start, self.end = self._get_endpoints(hsp1, hsp2, avg_insert_len)
        else:
            self.start, self.end = self.hsp1.hit_start, self.hsp1.hit_end

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
            f"matched={self.matched})"
        )

    def __len__(self):
        return self.end - self.start + 1

    def get_genes(self, buffer=4000):
        start = self.start - buffer
        end = self.end + buffer

        return [shift_feature(gene, start) for gene in self._genome[start:end].features]

    def _get_gene_coverage(self, gene):
        gene_start = gene.location.start
        gene_end = gene.location.end

        if gene_start >= self.start and gene_end <= self.end:
            return 1.0
        elif gene_end < self.start or gene_start > self.end:
            return 0.0
        elif gene_start < self.start and gene_end > self.end:  # TODO: check this
            return (self.end - self.start + 1) / (gene_end - gene_start + 1)
        elif gene_start < self.start and self.start <= gene_end <= self.end:
            return (gene_end - self.start + 1) / (gene_end - gene_start + 1)
        elif gene_end > self.end and self.start <= gene_start <= self.end:
            return (self.end - gene_start + 1) / (gene_end - gene_start + 1)
        else:
            raise RuntimeError(
                f"Unknown gene case: "
                f"{gene_start=}, {gene_end=}, "
                f"insert_start={self.start}, insert_end={self.end}"
            )

    def to_dataframe(self, buffer=4000):
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

    def get_graphic_records(
        self, buffer=4000, col1="#ebf3ed", col2="#2e8b57", feature_types=None, cmap=None
    ):
        genes = self.get_genes(buffer)
        if feature_types is None:
            feature_types = {x.type for x in genes}

        # Create a new graphic object for query sequence
        features = [
            GraphicFeature(
                start=self.start,
                end=self.end,
                strand=self.strand,
                color=col2,
                label=self.seq_id,
            )
        ]

        # Plot the query sequence on the upper axes
        record_seq = GraphicRecord(
            first_index=self.start - buffer,
            sequence_length=self.end - self.start + 2 * buffer,
            features=features,
        )

        # Create graphic objects for all the genes and CDSes using
        # dna_features_viewer.BiopythonTranslator()

        if cmap is not None:
            for i, gene in enumerate(genes):
                genes[i].qualifiers["color"] = cmap(self._get_gene_coverage(gene))

        conv = BiopythonTranslator()
        conv.default_feature_color = col1
        features = [conv.translate_feature(x) for x in genes if x.type in feature_types]

        # Plot the genes and CDSes in the region of the mapped sequence
        record_hits = GraphicRecord(
            first_index=self.start - buffer,
            sequence_length=self.end - self.start + 2 * buffer,
            features=features,
        )

        return record_seq, record_hits

    def plot(
        self,
        buffer=4000,
        col1="#ebf3ed",
        col2="#2e8b57",
        feature_types=None,
        colorbar=False,
        axs=None,
        backend="matplotlib",
        **kwargs,
    ):

        cmap = None
        if colorbar:
            cmap = LinearSegmentedColormap.from_list("custom", [col1, col2])

        seqs, hits = self.get_graphic_records(buffer, col1, col2, feature_types, cmap)
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
            _ = seqs.plot(ax=axs[0])
            _ = hits.plot(ax=axs[1])

            if colorbar:
                ax_cb = fig.add_axes([0.9, 0.1, 0.02, 0.8])
                fig.subplots_adjust(left=0.1, right=0.85)
                ColorbarBase(
                    ax_cb, cmap=cmap, orientation="vertical", label="gene coverage"
                )

            fig_axvline(axs, self.start)
            fig_axvline(axs, self.end)

            return axs
        else:
            raise NotImplementedError

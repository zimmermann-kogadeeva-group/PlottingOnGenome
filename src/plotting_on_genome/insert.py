import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import (
    BiopythonTranslator,
    GraphicFeature,
    GraphicRecord,
)

from .helper import shift_feature


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
        self, seq_id, seq1, hsp1, seq2=None, hsp2=None, genome=None, avg_insert_len=4000
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
        self.matched = True if hsp2 is not None else False

        self.start, self.end = self._get_endpoints(hsp1, hsp2, avg_insert_len)

        self._genome = genome

        self.cov1 = len(hsp1.query.seq) / len(seq1)
        self.cov2 = len(hsp2.query.seq) / len(seq2) if hsp2 is not None else None

        self.coverage = min(x for x in (self.cov1, self.cov2) if x is not None)

    def __len__(self):
        return self.end - self.start

    def get_genes(self, buffer=4000):
        start = self.start - buffer
        end = self.end + buffer

        return [shift_feature(gene, start) for gene in self._genome[start:end].features]

    def to_dataframe(self, buffer=4000):
        return pd.DataFrame(
            [
                {
                    "start": gene.location.start,
                    "end": gene.location.end,
                    "strand": gene.location.strand,
                    "type": gene.type,
                    **gene.qualifiers,
                }
                for gene in self.get_genes(buffer)
            ]
        )

    def _get_graphic_records_insert(self, buffer, col1, col2, feature_types=None):
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
        col1="#8DDEF7",
        col2="#CFFCCC",
        feature_types=None,
        figsize=None,
        axs=None,
        backend="matplotlib",
    ):

        seqs, hits = self._get_graphic_records_insert(buffer, col1, col2, feature_types)

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
        else:
            raise NotImplementedError

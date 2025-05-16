#!/usr/bin/env python3

import json
import re
from itertools import accumulate, product
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from dna_features_viewer import (
    CircularGraphicRecord,
    GraphicFeature,
    GraphicRecord,
)
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap

from .helper import get_genome, get_genome_file, run_blast, shift_feature
from .insert import Insert, fig_axvline


def _get_length(fwd, rev):
    """Helper function to calculate thee length of a potential insert before
    creating instantiating Insert class"""
    if fwd.hit_strand == +1:
        return rev.hit_end - fwd.hit_start
    else:
        return fwd.hit_end - rev.hit_start


def _get_contig_label(contig, mapped_ids, show_labels=True):
    if show_labels:
        if mapped_ids and contig.id in mapped_ids:
            return contig.id
        else:
            return None
    else:
        return None


class Mapping(object):

    def _get_paired_inserts(self, seq_id):
        # get all relevant fwd and rev hits from blast
        fwds, revs = [], []
        if seq_id + self.fwd_suf in self.blast_results:
            fwds = self.blast_results[seq_id + self.fwd_suf].hsps
        if seq_id + self.rev_suf in self.blast_results:
            revs = self.blast_results[seq_id + self.rev_suf].hsps

        # match appropriate hits
        matched_idxs = []
        unmatched_fwd_idxs = set(range(len(fwds)))
        unmatched_rev_idxs = set(range(len(revs)))

        for (i, fwd), (j, rev) in product(enumerate(fwds), enumerate(revs)):
            if (
                fwd.hit_id == rev.hit_id
                and fwd.hit_strand == -rev.hit_strand
                and 0 < _get_length(fwd, rev) < self.max_insert_len
            ):
                matched_idxs.append([i, j])
                unmatched_fwd_idxs.discard(i)
                unmatched_rev_idxs.discard(j)

        matched = [
            Insert(
                seq_id,
                idx,
                self.seqs[fwds[i].query_id],
                fwds[i],
                self.seqs[revs[j].query_id],
                revs[j],
                genome=self.genome[fwds[i].hit_id],
            )
            for idx, (i, j) in enumerate(matched_idxs)
        ]

        unmatched_fwd = [
            Insert(
                seq_id,
                idx + len(matched),
                self.seqs[fwds[i].query_id],
                fwds[i],
                genome=self.genome[fwds[i].hit_id],
                avg_insert_len=self.avg_insert_len,
            )
            for idx, i in enumerate(unmatched_fwd_idxs)
        ]

        unmatched_rev = [
            Insert(
                seq_id,
                idx + len(matched) + len(unmatched_fwd),
                self.seqs[revs[i].query_id],
                revs[i],
                genome=self.genome[revs[i].hit_id],
                avg_insert_len=self.avg_insert_len,
            )
            for idx, i in enumerate(unmatched_rev_idxs)
        ]

        return matched + unmatched_fwd + unmatched_rev

    def _get_single_inserts(self, seq_id):
        inserts = []
        if seq_id in self.blast_results:
            inserts = [
                Insert(
                    seq_id,
                    idx,
                    self.seqs[x.query_id],
                    x,
                    genome=self.genome[x.hit_id],
                    paired=False,
                )
                for idx, x in enumerate(self.blast_results[seq_id].hsps)
            ]
        return inserts

    def __init__(
        self,
        seq_file,
        work_dir,
        *,
        genome_file=None,
        search_term=None,
        email=None,
        retmax=200,
        fwd_suffix=None,
        rev_suffix=None,
        blast_clean=True,
        max_insert_len=10000,
        avg_insert_len=4000,
        **kwargs,
    ):
        # Populate obj attributes
        self.max_insert_len = max_insert_len
        self.avg_insert_len = avg_insert_len
        # If fwd and rev suffixes are None, then inserts are not paired
        self.fwd_suf = fwd_suffix
        self.rev_suf = rev_suffix

        # Make sure that specified work
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(exist_ok=True, parents=True)

        # Get seqs and remove primers if need be and redefine self.seq_file
        self.seq_file = Path(seq_file)
        self.seqs = SeqIO.to_dict(SeqIO.parse(self.seq_file, "fasta"))

        # Check at genome_file or search_term is specified
        if search_term is not None:
            self.genome_file = get_genome_file(self.work_dir, search_term, retmax)
            self.genome_fasta = self.genome_file.with_suffix(".fasta")
        elif genome_file is not None:
            self.genome_file = Path(genome_file)
            self.genome_fasta = self.work_dir / (self.genome_file.stem + ".fasta")
        else:
            raise ValueError("Either genome_file or search_term needs to be given")

        self.genome = get_genome(
            self.genome_file, self.genome_fasta, search_term, email, retmax
        )
        self.total_genome_length = sum(
            [len(x) for x in self.genome.values() if x.seq.defined]
        )

        # Save the input parameters
        with open(self.work_dir / "parameters.json", "w") as fh:
            json.dump(
                {
                    "seq_file": str(self.seq_file),
                    "search_term": search_term,
                    "genome_file": str(genome_file),
                    "email": email,
                    "work_dir": str(self.work_dir),
                    "retmax": retmax,
                    "fwd_suffix": self.fwd_suf,
                    "rev_suffix": self.rev_suf,
                },
                fh,
            )

        # Run BLAST and use xml format to save blast output
        self.blast_results = run_blast(
            self.seq_file,
            self.genome_fasta,
            self.work_dir / "blast_output.xml",
        )

        if blast_clean:
            self.blast_results = {
                k: v for k, v in self.blast_results.items() if len(v.hits) > 0
            }

        self.seqs = SeqIO.to_dict(SeqIO.parse(self.seq_file, "fasta"))

        if self.fwd_suf is not None and self.rev_suf is not None:
            self.seq_ids = sorted(
                {
                    re.sub(f"{self.fwd_suf}$|{self.rev_suf}$", "", x)
                    for x in self.seqs.keys()
                }
            )

            self._all_inserts = {
                seq_id: self._get_paired_inserts(seq_id) for seq_id in self.seq_ids
            }
        else:
            self.seq_ids = sorted(set(self.seqs.keys()))

            self._all_inserts = {
                seq_id: self._get_single_inserts(seq_id) for seq_id in self.seq_ids
            }

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._all_inserts[self.seq_ids[key]]

        elif isinstance(key, str):
            return self._all_inserts[key]

        elif isinstance(key, (tuple, list)):
            if all([isinstance(x, (str, int)) for x in key]):
                return [insert for k in key for insert in self.__getitem__(k)]
            elif all([isinstance(x, (tuple, list)) for x in key]):
                return [self._all_inserts[seq_id][ins_idx] for seq_id, ins_idx in key]
            else:
                raise ValueError(
                    "If key is tuple it is supposed to be list "
                    "of seq_ids or list of (seq_id, insert_idx) tuples"
                )

        elif isinstance(key, slice):
            return self.__getitem__(
                [ii for ii in range(*key.indices(len(self.seq_ids)))]
            )

        else:
            raise TypeError(f"Invalid argument type: {type(key)}")

    def __len__(self):
        return len(self.seq_ids)

    def get_genes(self, start, end, hit_id, buffer=4000):
        start_ = start - buffer
        end_ = end + buffer

        return [
            shift_feature(gene, start_)
            for gene in self.genome[hit_id][start_:end_].features
        ]

    def get(
        self,
        seq_id_or_idx=None,
        insert_type="both",
        filter_threshold=None,
    ):
        assert insert_type in ("matched", "unmatched", "both"), f"{insert_type=}"
        if seq_id_or_idx is None:
            seq_id_or_idx = self.seq_ids

        inserts = self.__getitem__(seq_id_or_idx)

        # Apply the coverage filter
        if filter_threshold is not None:
            if not 0 <= filter_threshold <= 1:
                raise ValueError("Filter value needs to be between 0 and 1")

            inserts = [x for x in inserts if x.coverage > filter_threshold]

        if insert_type == "matched":
            inserts = [x for x in inserts if x.matched]
        elif insert_type == "unmatched":
            inserts = [x for x in inserts if not x.matched]

        return inserts

    def get_insert_ids(
        self, seq_id_or_idx=None, insert_type="both", filter_threshold=None, **kwargs
    ):
        return [
            (x.seq_id, x.idx)
            for x in self.get(
                seq_id_or_idx=seq_id_or_idx,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
            )
        ]

    def get_seq_ids(self, insert_type="both", filter_threshold=None, **kwargs):
        return [
            x.seq_id
            for x in self.get(
                insert_type=insert_type, filter_threshold=filter_threshold
            )
        ]

    def to_dataframe(self, insert_type="both", filter_threshold=None):
        return pd.DataFrame(
            [
                (
                    x.seq_id,
                    x.hit_id,
                    x.idx,
                    x.start,
                    x.end,
                    x.strand,
                    len(x),
                    x.coverage,
                    x.matched,
                )
                for x in self.get(
                    insert_type=insert_type, filter_threshold=filter_threshold
                )
            ],
            columns=(
                "seq_id",
                "hit_id",
                "insert_idx",
                "insert_start",
                "insert_end",
                "insert_strand",
                "insert_length",
                "insert_coverage",
                "insert_matched",
            ),
        )

    def genes_to_dataframe(
        self,
        seq_id_or_idx=None,
        insert_type="both",
        filter_threshold=None,
        buffer=4000,
    ):
        inserts = self.get(seq_id_or_idx, insert_type, filter_threshold)

        df_genes = pd.DataFrame(
            [], columns=("start", "end", "strand", "type", "coverage", "locus_tag")
        )
        if len(inserts):
            df_genes = pd.concat(
                [
                    insert.to_dataframe(buffer).assign(
                        insert_idx=insert.idx, seq_id=insert.seq_id
                    )
                    for insert in inserts
                ]
            ).reset_index(drop=True)

        return df_genes

    # TODO: split this into two: one for genome features and one for insert features
    def get_graphic_features(
        self,
        seq_id_or_idxs=None,
        insert_type="both",
        filter_threshold=None,
        show_labels=True,
        col1="#ebf3ed",
        col2="#2e8b57",
        **kwargs,
    ):

        inserts = self.get(seq_id_or_idxs, insert_type, filter_threshold)

        # Get just the sequences for each NCBI record and order them by size in
        # descending order. 'x.features[0]' to get the whole sequence for a
        # given NCBI record. Other features are specific genes, cfs, etc.
        db_seqs = sorted(
            [x for x in self.genome.values() if x.seq.defined],
            key=lambda x: len(x),
            reverse=True,
        )

        linecolor = "#000000"
        if "linecolor" in kwargs:
            linecolor = kwargs.pop("linecolor")

        # Get the shifts needed to plot all the NCBI records in a continuous line
        shifts = list(accumulate([0] + [len(x) for x in db_seqs]))

        # Get IDs of NCBI records that were mapped to. Used to check where to
        # add labels if option is set.
        mapped_ids = {x.hit_id for x in inserts}

        # Make plots of NCBI records and label only the ones that were mapped
        # to. Using BiopythonTranslator() didn't allow for control of labels,
        # hence we are just using GraphicFeature class
        genome = [
            (
                GraphicFeature(
                    start=shifts[i],
                    end=shifts[i + 1],
                    label=_get_contig_label(x, mapped_ids, show_labels),
                    color=col1,
                    linecolor="#000000",
                    **kwargs,
                )
            )
            for i, x in enumerate(db_seqs)
        ]

        # Get IDs of NCBI records in the order as in the figure. Used to make
        # sure locations are shifted correctly.
        ids = [x.id for x in db_seqs]

        # Add plots of the query sequences plotted on top of the plots of NCBI records
        hits = [
            GraphicFeature(
                start=insert.start + shifts[ids.index(insert.hit_id)] + 1,
                end=insert.end + shifts[ids.index(insert.hit_id)],
                strand=insert.strand,
                color=col2,
                label=None,
                linecolor=linecolor,
                **kwargs,
            )
            for insert in inserts
        ]

        return shifts[-1], genome + hits

    def plot(
        self,
        seq_id_or_idxs=None,
        insert_type="both",
        filter_threshold=None,
        show_labels=True,
        col1="#ebf3ed",
        col2="#2e8b57",
        ax=None,
        backend="matplotlib",
        **kwargs,
    ):
        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        seq_len, features = self.get_graphic_features(
            seq_id_or_idxs,
            insert_type,
            filter_threshold,
            show_labels,
            col1,
            col2,
            **kwargs,
        )

        rec = CircularGraphicRecord(sequence_length=seq_len, features=features)

        if backend == "matplotlib":
            if ax is None:
                fig, ax = plt.subplots(1, 1, **kwargs)

            _ = rec.plot(ax, annotate_inline=False)
        else:
            raise NotImplementedError

        return ax

    def get_clusters(self, insert_type="both", filter_threshold=None):
        inserts = sorted(
            self.get(insert_type=insert_type, filter_threshold=filter_threshold),
            key=lambda x: (x.start, x.end),
        )

        if len(inserts):
            clusters = []
            current_cluster = [inserts[0]]
            for i in range(1, len(inserts)):
                current_insert = inserts[i]
                last_insert = current_cluster[-1]

                # Check if the current interval overlaps with the last interval in the
                # cluster
                if current_insert.start <= last_insert.end:
                    current_cluster.append(current_insert)
                else:
                    # If no overlap, start a new cluster
                    clusters.append(current_cluster)
                    current_cluster = [current_insert]

            # Add the last cluster
            clusters.append(current_cluster)

            return [
                [(x.seq_id, x.idx) for x in cluster]
                for cluster in clusters
                if len(cluster) > 1
            ]
        else:
            return []

    def plot_inserts(
        self,
        seq_id_or_idx,
        insert_type="both",
        filter_threshold=None,
        buffer=4000,
        col1="#ebf3ed",
        col2="#2e8b57",
        feature_types=None,
        colorbar=False,
        axs=None,
        backend="matplotlib",
        **kwargs,
    ):

        inserts = self.get(seq_id_or_idx, insert_type, filter_threshold)
        cmap = None
        if colorbar:
            cmap = LinearSegmentedColormap.from_list("custom", [col1, col2])

        features = [x.get_graphic_feature(col2) for x in inserts]

        # TODO: uncomment and use get_genes from InsertsDict
        start = inserts[0].start  # min([x.start for x in inserts])
        end = inserts[0].end  # max([x.end for x in inserts])

        # Plot the query sequence on the upper axes
        rec_seqs = GraphicRecord(
            first_index=start - buffer,
            sequence_length=end - start + 2 * buffer,
            features=features,
        )
        rec_genes = inserts[0].get_genes_graphic_record(
            buffer, col1, feature_types, cmap
        )

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
            _ = rec_seqs.plot(ax=axs[0])
            _ = rec_genes.plot(ax=axs[1])

            if colorbar:
                ax_cb = fig.add_axes([0.9, 0.1, 0.02, 0.8])
                fig.subplots_adjust(left=0.1, right=0.85)
                ColorbarBase(
                    ax_cb, cmap=cmap, orientation="vertical", label="gene coverage"
                )

            fig_axvline(axs, inserts[0].start)
            fig_axvline(axs, inserts[-1].end)

            return axs
        else:
            raise NotImplementedError

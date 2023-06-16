#!/usr/bin/env python3

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import subprocess

from Bio import Entrez, SeqIO, SearchIO
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    BiopythonTranslator,
    CircularGraphicRecord,
)


def correct_hit_id(x):
    x.id = x.id.split("|")[1]
    return x


def shift_feature(feature, shift=0):
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def set_feature(feature, **kwargs):
    new_feature = deepcopy(feature)
    for k, v in kwargs.items():
        if hasattr(new_feature, k):
            setattr(new_feature, k, v)
    return new_feature


def get_db(search_term, email, gbk_file=None, retmax=100):
    # TODO: maybe replace checking for existing file with caching in user's home
    if gbk_file is not None and Path(gbk_file).exists():
        # Read in the file
        data_gb = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))

    else:
        # Set the email address for NCBI queries
        Entrez.email = email

        # Get record ids from NCBI
        with Entrez.esearch(db="nucleotide", term=search_term, retmax=retmax) as handle:
            res = Entrez.read(handle)

        # Download these records from NCBI rettype="gbwithparts" is used to
        # download the sequences with record features
        with Entrez.efetch(
            db="nucleotide",
            id=",".join(res["IdList"]),
            rettype="gbwithparts",
            retmode="text",
        ) as handle:
            data_gb = SeqIO.to_dict(SeqIO.parse(handle, "gb"))

        # Save the file if it doesn't exist
        if gbk_file is not None and not Path(gbk_file).exists():
            path = Path(gbk_file)

            # Save in genbank format (stores more information about the NCBI records)
            SeqIO.write(data_gb.values(), path, "genbank")

            # Save in fasta format (only acceptable format for makeblastdb)
            SeqIO.write(
                [x for x in data_gb.values() if x.seq.defined],
                path.with_suffix(".fasta"),
                "fasta",
            )

    return data_gb


def run_blast(seq_file, db_file, blast_output):
    # make blast database
    subprocess.run(
        f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # align input sequences with query using blastn with default options
    subprocess.run(
        f"blastn -query {seq_file} -db {db_file} -out {blast_output} -outfmt 5",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    return {
        x.id: x.hit_map(correct_hit_id)
        for x in SearchIO.parse(blast_output, "blast-xml")
    }


class Pipeline(object):
    def __init__(
        self,
        seq_file,
        search_term,
        email,
        work_dir,
        retmax=100,
        fwd_suffix="_F",
        rev_suffix="_R",
        blast_clean=True,
        **kwargs,
    ):
        self.seq_file = seq_file
        self.work_dir = Path(work_dir)
        self.search_term = search_term
        self.__strand_func__ = lambda seq_id: +1 if seq_id.endswith(fwd_suffix) else -1

        self._db_col = "#8DDEF7" if "db_col" not in kwargs else kwargs["db_col"]
        self._query_col = (
            "#CFFCCC" if "query_col" not in kwargs else kwargs["query_col"]
        )

        self.seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

        self.db = get_db(search_term, email, self.work_dir / "db.gbk", retmax)

        # Use xml format to save blast output
        self.blast_results = run_blast(
            seq_file, self.work_dir / "db.fasta", self.work_dir / "blast_output.xml"
        )

        if blast_clean:
            self.blast_results = {
                k: v for k, v in self.blast_results.items() if len(v.hits) > 0
            }

    def plot_hsp(self, seq_id, hsp_idx=0, buffer=4000, figsize=None, save_fmt=None):
        # Default values for figure size
        figsize = figsize or (10, 8)

        # Get the specified HSP object from dictionary of BLAST results
        hsp = self.blast_results[seq_id].hsps[hsp_idx]

        seq_mapped_start = hsp.hit_start - hsp.query_start
        seq_mapped_end = hsp.hit_end + (len(self.seqs[hsp.query_id]) - hsp.query_end)

        # Create a new figure and axes object
        fig, axs = plt.subplots(2, 1, figsize=figsize)

        # Create a new graphic object for query sequence
        features = [
            GraphicFeature(
                start=seq_mapped_start,
                end=seq_mapped_end,
                strand=self.__strand_func__(seq_id),
                color=self._query_col,
                label=hsp.query_id,
            )
        ]

        # Plot the query sequence on the upper axes
        record_seq = GraphicRecord(
            first_index=seq_mapped_start - buffer,
            sequence_length=seq_mapped_end - seq_mapped_start + 2 * buffer,
            features=features,
        )
        _ = record_seq.plot(ax=axs[0])

        # Create graphic objects for all the genes and CDSes using
        # dna_features_viewer.BiopythonTranslator()
        conv = BiopythonTranslator()
        conv.default_feature_color = self._db_col
        features = [
            conv.translate_feature(shift_feature(x, seq_mapped_start - buffer))
            for x in self.db[hsp.hit_id][
                seq_mapped_start - buffer : seq_mapped_end + buffer
            ].features
        ]

        # Plot the genes and CDSes in the region of the mapped sequence
        record_hits = GraphicRecord(
            first_index=seq_mapped_start - buffer,
            sequence_length=seq_mapped_end - seq_mapped_start + 2 * buffer,
            features=features,
        )
        _ = record_hits.plot(ax=axs[1])

        if save_fmt is not None:
            fig.savefig(self.work_dir / f"{seq_id}_hit{hsp_idx}.{save_fmt}")

        return fig

    def plot_all_hsp(self, seq_id, buffer=4000, figsize=None, save_fmt=None):
        return [
            self.plot_hsp(seq_id, i, buffer, figsize, save_fmt)
            for i, x in enumerate(self.blast_results[seq_id].hsps)
        ]

    def plot_all_db_seqs(self, labels=True, figsize=None, save_fmt=None):
        figsize = figsize or (10, 30)
        fig, ax = plt.subplots(figsize=figsize)
        # Get just the sequences for each NCBI record and order them by size in
        # descending order. 'x.features[0]' to get the whole sequence for a
        # given NCBI record. Other features are specific genes, cfs, etc.
        db_seqs = sorted(
            [
                set_feature(x.features[0], id=x.id)
                for x in self.db.values()
                if x.seq.defined
            ],
            key=lambda x: len(x),
            reverse=True,
        )

        # Get the shifts needed to plot all the NCBI records in a continuous line
        shifts = np.cumsum([0] + [len(x) for x in db_seqs])

        # Get IDs of NCBI records that were mapped to
        mapped_ids = set([x.id for l in self.blast_results.values() for x in l.hits])

        # Make plots of NCBI records and label only the ones that were mapped
        # to. Using BiopythonTranslator() didn't allow for control of labels,
        # hence we are just using GraphicFeature class
        features = [
            GraphicFeature(
                start=x.location.start + shifts[i],
                end=x.location.end + shifts[i],
                strand=x.strand,
                label=x.id if labels else None,
                color=self._db_col,
            )
            if x.id in mapped_ids
            else GraphicFeature(
                start=x.location.start + shifts[i],
                end=x.location.end + shifts[i],
                strand=x.strand,
                color=self._db_col,
            )
            for i, x in enumerate(db_seqs)
        ]

        # Get IDs of NCBI records in the order as in the figure. Used to make
        # sure locations are shifted correctly.
        ids = [x.id for x in db_seqs]

        # Add plots of the query sequences plotted on top of the plots of NCBI records
        hits = [
            GraphicFeature(
                start=x.hit_start - x.query_start + shifts[ids.index(x.hit_id)],
                end=x.hit_end
                + (len(self.seqs[x.query_id]) - x.query_end)
                + shifts[ids.index(x.hit_id)],
                strand=self.__strand_func__(x.query_id),
                color=self._query_col,
                label=x.query_id if labels else None,
            )
            for l in self.blast_results.values()
            for x in l.hsps
        ]

        rec = CircularGraphicRecord(
            sequence_length=shifts[-1], features=features + hits
        )
        _ = rec.plot(ax, annotate_inline=False)

        if save_fmt is not None:
            fig.savefig(self.work_dir / f"genome_plot.{save_fmt}")

        return fig

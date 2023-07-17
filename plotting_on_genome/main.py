#!/usr/bin/env python3

from copy import deepcopy
import hashlib
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re
import subprocess

from Bio import Entrez, SeqIO, SearchIO
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    BiopythonTranslator,
    CircularGraphicRecord,
)


def correct_hit_id(x):
    """Helper function to get rid red| or emb| added by BLAST to contig ID"""
    x.id = x.id.split("|")[1]
    return x


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def set_feature(feature, **kwargs):
    """Helper function to set an attribute in Biopython feature"""
    new_feature = deepcopy(feature)
    for k, v in kwargs.items():
        if hasattr(new_feature, k):
            setattr(new_feature, k, v)
    return new_feature


def get_length(fwd, rev):
    """Helper function to calculate thee length of a potential insert before
    creating instantiating Insert class"""
    if fwd.hit_strand == +1:
        return rev.hit_end - fwd.hit_start
    else:
        return fwd.hit_end - rev.hit_start


def get_db(search_term, email, retmax=100, out_prefix=None):
    gbk_file = None
    if out_prefix is not None:
        gbk_file = Path(out_prefix).with_suffix(".gbk")
        fasta_file = Path(out_prefix).with_suffix(".fasta")

    if gbk_file is not None and Path(gbk_file).exists():
        # Read in the file
        return SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))

    # Otherwise retrieve the records from NCBI
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
            # Save in genbank format (stores more information about the NCBI records)
            SeqIO.write(data_gb.values(), gbk_file, "genbank")

            # Save in fasta format (only acceptable format for makeblastdb)
            SeqIO.write(
                [x for x in data_gb.values() if x.seq.defined],
                fasta_file,
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

    # Return a dictionary correcting the hit IDs which get an unnecessary
    # prefix from BLAST
    return {
        x.id: x.hit_map(correct_hit_id)
        for x in SearchIO.parse(blast_output, "blast-xml")
    }


class Insert(object):
    def __init__(self, hsp1, hsp2=None, avg_insert_length=4000):
        # TODO: maybe switch to using a list of hsps instead of hsp1 and hsp2
        self.hsp1 = hsp1
        self.hsp2 = hsp2
        self.strand = self.hsp1.hit_strand
        self.hit_id = self.hsp1.hit_id
        self.query_id = self.hsp1.query_id

        if hsp2 is not None:
            if self.hsp1.hit_strand == 1:
                self.start = self.hsp1.hit_start
                self.end = self.hsp2.hit_end
            else:
                self.start = self.hsp2.hit_start
                self.end = self.hsp1.hit_end
        else:
            self.start = self.hsp1.hit_start
            self.end = self.hsp1.hit_start + avg_insert_length

    def __len__(self):
        return self.end - self.start


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
        self.__fwd_suf__ = fwd_suffix
        self.__rev_suf__ = rev_suffix

        self.work_dir.mkdir(exist_ok=True)

        self._db_col = "#8DDEF7" if "db_col" not in kwargs else kwargs["db_col"]
        self._query_col = (
            "#CFFCCC" if "query_col" not in kwargs else kwargs["query_col"]
        )

        self.seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

        self.seq_ids = list(
            {
                re.sub(f"{self.__fwd_suf__}$|{self.__rev_suf__}$", "", x)
                for x in self.seqs.keys()
            }
        )

        # Hash the search term to use as filename in cache dir
        search_hashed = hashlib.sha1(search_term.encode()).hexdigest()
        db_path = self.work_dir / f"db_{search_hashed}"
        self.db = get_db(search_term, email, retmax, db_path)

        # Run BLAST and use xml format to save blast output
        self.blast_results = run_blast(
            seq_file, db_path.with_suffix(".fasta"), self.work_dir / "blast_output.xml"
        )

        if blast_clean:
            self.blast_results = {
                k: v for k, v in self.blast_results.items() if len(v.hits) > 0
            }

    def get_inserts(
        self, seq_id, output="both", filter_threshold=None, insert_max_len=10000
    ):
        matched = []

        fwds, revs = [], []
        if seq_id + self.__fwd_suf__ in self.blast_results:
            fwds = self.blast_results[seq_id + self.__fwd_suf__].hsps
        if seq_id + self.__rev_suf__ in self.blast_results:
            revs = self.blast_results[seq_id + self.__rev_suf__].hsps

        if filter_threshold is not None:
            assert 0 <= filter_threshold <= 1
            fwds = [
                x
                for x in fwds
                if len(x.query.seq) / len(self.seqs[x.query_id]) > filter_threshold
            ]
            revs = [
                x
                for x in revs
                if len(x.query.seq) / len(self.seqs[x.query_id]) > filter_threshold
            ]

        unmatched = set(fwds).union(revs)

        for fwd, rev in product(fwds, revs):
            if (
                fwd.hit_id == rev.hit_id
                and fwd.hit_strand == -rev.hit_strand
                and 0 < get_length(fwd, rev) < insert_max_len
            ):
                matched.append([fwd, rev])
                unmatched -= {fwd, rev}

        matched = [Insert(*x) for x in matched]
        unmatched = [Insert(x) for x in unmatched]

        if output == "matched":
            return matched
        elif output == "unmatched":
            return unmatched
        else:
            return matched + unmatched

    def plot_insert(self, insert, buffer=4000, figsize=None, axs=None):
        # Default values for figure size and create the figure
        if axs is None:
            figsize = figsize or (10, 8)
            fig, axs = plt.subplots(2, 1, figsize=figsize)
        assert len(axs) == 2

        # Get the sequence id
        seq_id = re.sub(f"{self.__fwd_suf__}|{self.__rev_suf__}$", "", insert.query_id)

        # Create a new graphic object for query sequence
        features = [
            GraphicFeature(
                start=insert.start + 1,  # +1 due to python's 0-indexing
                end=insert.end,
                strand=insert.strand,
                color=self._query_col,
                label=seq_id,
            )
        ]

        # Plot the query sequence on the upper axes
        record_seq = GraphicRecord(
            first_index=insert.start - buffer + 1,  # +1 due to python's 0-indexing
            sequence_length=insert.end - insert.start + 2 * buffer,
            features=features,
        )
        _ = record_seq.plot(ax=axs[0])

        # Create graphic objects for all the genes and CDSes using
        # dna_features_viewer.BiopythonTranslator()
        conv = BiopythonTranslator()
        conv.default_feature_color = self._db_col
        features = [
            conv.translate_feature(shift_feature(x, insert.start - buffer))
            for x in self.db[insert.hit_id][
                insert.start - buffer : insert.end + buffer
            ].features
        ]

        # Plot the genes and CDSes in the region of the mapped sequence
        record_hits = GraphicRecord(
            first_index=insert.start - buffer + 1,  # +1 due to python's 0-indexing
            sequence_length=insert.end - insert.start + 2 * buffer,
            features=features,
        )
        _ = record_hits.plot(ax=axs[1])

        return axs

    def plot_all_inserts(
        self,
        seq_id,
        output="both",
        filter_threshold=None,
        insert_max_len=10000,
        buffer=4000,
        figsize=None,
    ):
        return [
            self.plot_insert(x, buffer, figsize)
            for i, x in enumerate(
                self.get_inserts(
                    seq_id,
                    output=output,
                    filter_threshold=filter_threshold,
                    insert_max_len=insert_max_len,
                )
            )
        ]

    def plot_all_db_seqs(
        self,
        output="both",
        filter_threshold=None,
        insert_max_len=10000,
        labels=True,
        figsize=None,
        ax=None,
    ):
        if ax is None:
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

        # Get IDs of NCBI records that were mapped to. Used to check where to
        # add labels if option is set.
        mapped_ids = set([x.id for l in self.blast_results.values() for x in l.hits])

        # Make plots of NCBI records and label only the ones that were mapped
        # to. Using BiopythonTranslator() didn't allow for control of labels,
        # hence we are just using GraphicFeature class
        features = [
            GraphicFeature(
                start=x.location.start + shifts[i] + 1,  # +1 due to python's 0-indexing
                end=x.location.end + shifts[i],
                strand=x.strand,
                label=x.id if labels else None,
                color=self._db_col,
            )
            if x.id in mapped_ids
            else GraphicFeature(
                start=x.location.start + shifts[i] + 1,  # +1 due to python's 0-indexing
                end=x.location.end + shifts[i],
                strand=x.strand,
                color=self._db_col,
            )
            for i, x in enumerate(db_seqs)
        ]

        # Get IDs of NCBI records in the order as in the figure. Used to make
        # sure locations are shifted correctly.
        ids = [x.id for x in db_seqs]
        seq_ids = {
            re.sub(f"{self.__fwd_suf__}|{self.__rev_suf__}$", "", x)
            for x in self.blast_results.keys()
        }

        # Add plots of the query sequences plotted on top of the plots of NCBI records
        hits = [
            GraphicFeature(
                start=insert.start + shifts[ids.index(insert.hit_id)] + 1,
                end=insert.end + shifts[ids.index(insert.hit_id)],
                strand=insert.strand,
                color=self._query_col,
                label=None,
            )
            for seq_id in seq_ids
            for insert in self.get_inserts(
                seq_id,
                output=output,
                filter_threshold=filter_threshold,
                insert_max_len=insert_max_len,
            )
        ]

        rec = CircularGraphicRecord(
            sequence_length=shifts[-1], features=features + hits
        )

        _ = rec.plot(ax, annotate_inline=False)

        return ax

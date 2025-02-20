#!/usr/bin/env python3

import hashlib
import json
import re
from itertools import accumulate, product
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from dna_features_viewer import (
    CircularGraphicRecord,
    GraphicFeature,
)

from .helper import download_genome, run_blast, shift_feature
from .insert import Insert


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


class InsertsDict(object):

    def _get_genome_file(self, work_dir, search_term, retmax):
        # Hash the search term to use as filename in cache dir
        search_hashed = hashlib.sha1((search_term + str(retmax)).encode()).hexdigest()
        return work_dir / f"db_{search_hashed}.gbk"

    def _get_genome(
        self, genome_file, genome_fasta, search_term=None, email=None, retmax=None
    ):
        # Get genome
        if search_term is not None:
            if email is None:
                raise RuntimeError("Email is required for NCBI API")

            genome = download_genome(search_term, email, retmax, genome_file)
        else:
            if genome_file.suffix == ".gff":
                genome = SeqIO.to_dict(GFF.parse(genome_file))
            elif genome_file.suffix == ".gbk":
                genome = SeqIO.to_dict(SeqIO.parse(genome_file, "genbank"))
            else:
                raise RuntimeError(
                    "Wrong format expected either `.gbk` or `.gff` file."
                )

        # Save in fasta format (only acceptable format for makeblastdb)
        SeqIO.write(
            [x for x in genome.values() if x.seq.defined],
            genome_fasta,
            "fasta",
        )

        return genome

    def _get_paired_inserts(self, seq_id):
        # get all relevant fwd and rev hits from blast
        fwds, revs = [], []
        if seq_id + self._fwd_suf in self._blast_results:
            fwds = self._blast_results[seq_id + self._fwd_suf].hsps
        if seq_id + self._rev_suf in self._blast_results:
            revs = self._blast_results[seq_id + self._rev_suf].hsps

        # match appropriate hits
        matched_idxs = []
        unmatched_fwd_idxs = set(range(len(fwds)))
        unmatched_rev_idxs = set(range(len(revs)))

        for (i, fwd), (j, rev) in product(enumerate(fwds), enumerate(revs)):
            if (
                fwd.hit_id == rev.hit_id
                and fwd.hit_strand == -rev.hit_strand
                and 0 < _get_length(fwd, rev) < self._max_insert_len
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
                avg_insert_len=self._avg_insert_len,
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
                avg_insert_len=self._avg_insert_len,
            )
            for idx, i in enumerate(unmatched_rev_idxs)
        ]

        return matched + unmatched_fwd + unmatched_rev

    def _get_single_inserts(self, seq_id):
        return [
            Insert(
                seq_id,
                idx,
                self.seqs[x.query_id],
                x,
                genome=self.genome[x.hit_id],
                avg_insert_len=0,  # TODO: check that it makes sense
            )
            for idx, x in enumerate(self._blast_results[seq_id].hsps)
        ]

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
        self.seq_file = seq_file
        self.work_dir = Path(work_dir)
        self._max_insert_len = max_insert_len
        self._avg_insert_len = avg_insert_len
        # If fwd and rev suffixes are None, then inserts are not paired
        self._fwd_suf = fwd_suffix
        self._rev_suf = rev_suffix

        # Make sure that specified work
        self.work_dir.mkdir(exist_ok=True, parents=True)

        # Check at genome_file or search_term is specified
        if search_term is not None:
            genome_file = self._get_genome_file(self.work_dir, search_term, retmax)
            genome_fasta = genome_file.with_suffix(".fasta")
        elif genome_file is not None:
            genome_file = Path(genome_file)
            genome_fasta = self.work_dir / (genome_file.stem + ".fasta")
        else:
            raise ValueError("Either genome_file or search_term needs to be given")

        self._genome = self._get_genome(
            genome_file, genome_fasta, search_term, email, retmax
        )

        # Save the input parameters
        with open(self.work_dir / "parameters.json", "w") as fh:
            json.dump(
                {
                    "seq_file": self.seq_file,
                    "search_term": search_term,
                    "genome_file": str(genome_file),
                    "email": email,
                    "work_dir": str(self.work_dir),
                    "retmax": retmax,
                    "fwd_suffix": self._fwd_suf,
                    "rev_suffix": self._rev_suf,
                },
                fh,
            )

        # Run BLAST and use xml format to save blast output
        self._blast_results = run_blast(
            seq_file,
            genome_fasta,
            self.work_dir / "blast_output.xml",
        )

        if blast_clean:
            self._blast_results = {
                k: v for k, v in self._blast_results.items() if len(v.hits) > 0
            }

        self._seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

        if self._fwd_suf is not None and self._rev_suf is not None:
            self._seq_ids = tuple(
                {
                    re.sub(f"{self._fwd_suf}$|{self._rev_suf}$", "", x)
                    for x in self.seqs.keys()
                }
            )

            self._all_inserts = {
                seq_id: self._get_paired_inserts(seq_id) for seq_id in self._seq_ids
            }
        else:
            self._seq_ids = set(self.seqs.keys())

            self._all_inserts = {
                seq_id: self._get_single_inserts(seq_id) for seq_id in self._seq_ids
            }

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._all_inserts[self._seq_ids[key]]

        elif isinstance(key, str):
            return self._all_inserts[key]

        elif isinstance(key, (tuple, list)):
            return [insert for k in key for insert in self.__getitem__(k)]

        elif isinstance(key, slice):
            return self.__getitem__(
                [ii for ii in range(*key.indices(len(self.seq_ids)))]
            )
        else:
            raise TypeError(f"Invalid argument type: {type(key)}")

    def __len__(self):
        return len(self._seq_ids)

    @property
    def seqs(self):
        return self._seqs

    @property
    def seq_ids(self):
        return self._seq_ids

    @property
    def genome(self):
        return self._genome

    @property
    def blast_results(self):
        return self._blast_results

    def get_genes(self, start, end, hit_id, buffer=4000):
        start_ = start - buffer
        end_ = end + buffer

        return [
            shift_feature(gene, start_)
            for gene in self._genome[hit_id][start_:end_].features
        ]

    def get(
        self,
        seq_id_or_idx=None,
        insert_types="both",
        filter_threshold=None,
    ):
        assert insert_types in ("matched", "unmatched", "both")
        seq_id_or_idx = seq_id_or_idx or self.seq_ids

        inserts = self.__getitem__(seq_id_or_idx)
        # Apply the coverage filter
        if filter_threshold is not None:
            if not 0 <= filter_threshold <= 1:
                raise ValueError("Filter value needs to be between 0 and 1")

            inserts = [x for x in inserts if x.coverage > filter_threshold]

        if insert_types == "matched":
            inserts = [x for x in inserts if x.matched]
        elif insert_types == "unmatched":
            inserts = [x for x in inserts if not x.matched]

        return inserts

    def to_dataframe(self, insert_types="both", filter_threshold=None):
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
                )
                for x in self.get(
                    insert_types=insert_types, filter_threshold=filter_threshold
                )
            ],
            columns=(
                "seq_id",
                "NCBI_accession_number",
                "insert_idx",
                "insert_start",
                "insert_end",
                "insert_strand",
                "insert_length",
                "insert_coverage",
            ),
        )

    def genes_to_dataframe(
        self,
        seq_id_or_idx=None,
        insert_types="both",
        filter_threshold=None,
        buffer=4000,
    ):
        inserts = self.get(seq_id_or_idx, insert_types, filter_threshold)

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

    def _get_graphic_records_genome(self, inserts, show_labels, col1, col2):

        # Get just the sequences for each NCBI record and order them by size in
        # descending order. 'x.features[0]' to get the whole sequence for a
        # given NCBI record. Other features are specific genes, cfs, etc.
        db_seqs = sorted(
            [x for x in self._genome.values() if x.seq.defined],
            key=lambda x: len(x),
            reverse=True,
        )

        # Get the shifts needed to plot all the NCBI records in a continuous line
        shifts = list(accumulate([0] + [len(x) for x in db_seqs]))

        # Get IDs of NCBI records that were mapped to. Used to check where to
        # add labels if option is set.
        mapped_ids = {x.hit_id for x in inserts}

        # Make plots of NCBI records and label only the ones that were mapped
        # to. Using BiopythonTranslator() didn't allow for control of labels,
        # hence we are just using GraphicFeature class
        features = [
            (
                GraphicFeature(
                    start=shifts[i],
                    end=shifts[i + 1],
                    label=_get_contig_label(x, mapped_ids, show_labels),
                    color=col1,
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
            )
            for insert in inserts
        ]

        rec = CircularGraphicRecord(
            sequence_length=shifts[-1], features=features + hits
        )

        return rec

    def plot(
        self,
        inserts=None,
        show_labels=True,
        col1="#8DDEF7",
        col2="#CFFCCC",
        ax=None,
        backend="matplotlib",
        **kwargs,
    ):
        if "figsize" not in kwargs:
            kwargs["figsize"] = (10, 8)

        inserts = inserts or []
        rec = self._get_graphic_records_genome(inserts, show_labels, col1, col2)

        if backend == "matplotlib":
            if ax is None:
                fig, ax = plt.subplots(1, 1, **kwargs)

            _ = rec.plot(ax, annotate_inline=False)
        else:
            raise NotImplementedError

        return ax

#!/usr/bin/env python3

import hashlib
import json
import re
import subprocess
from itertools import product
from pathlib import Path

import pandas as pd
from BCBio import GFF
from Bio import Entrez, SearchIO, SeqIO


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def correct_hit_id(x):
    """Helper function to get rid red| or emb| added by BLAST to contig ID"""
    if "|" in x.id:
        x.id = x.id.split("|")[1]
    return x


def get_length(fwd, rev):
    """Helper function to calculate thee length of a potential insert before
    creating instantiating Insert class"""
    if fwd.hit_strand == +1:
        return rev.hit_end - fwd.hit_start
    else:
        return fwd.hit_end - rev.hit_start


def download_genome(search_term, email, retmax=100, gbk_file=None):
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
    def __init__(
        self, hsp1, hsp2=None, avg_insert_length=4000, seq_id=None, seq_len=None
    ):
        # TODO: maybe switch to using a list of hsps instead of hsp1 and hsp2
        self.hsp1 = hsp1
        self.hsp2 = hsp2
        self.strand = self.hsp1.hit_strand
        self.hit_id = self.hsp1.hit_id
        self.query_id = self.hsp1.query_id
        self.seq_id = seq_id
        self.seq_len = seq_len
        self.coverage = None

        # TODO: check seq_len
        if hsp2 is not None:
            if self.hsp1.hit_strand == 1:
                self.start = self.hsp1.hit_start
                self.end = self.hsp2.hit_end
            else:
                self.start = self.hsp2.hit_start
                self.end = self.hsp1.hit_end

            self.coverage = [
                len(self.hsp1.query.seq) / seq_len[0],
                len(self.hsp2.query.seq) / seq_len[1],
            ]
        else:
            self.start = self.hsp1.hit_start
            self.end = self.hsp1.hit_start + avg_insert_length
            self.coverage = [len(self.hsp1.query.seq) / seq_len[0]]

    def __len__(self):
        return self.end - self.start


class Pipeline(object):
    def __init__(
        self,
        seq_file,
        work_dir,
        genome_file=None,
        search_term=None,
        email=None,
        retmax=200,
        fwd_suffix="_F",
        rev_suffix="_R",
        blast_clean=True,
        **kwargs,
    ):
        # Check at genome_file or search_term is specified
        assert (
            genome_file is not None or search_term is not None
        ), "Either genome_file or search_term needs to be given"

        # Populate class attributes
        self.seq_file = seq_file
        self.work_dir = Path(work_dir)
        self._fwd_suf = fwd_suffix
        self._rev_suf = rev_suffix

        # Make sure that specified work
        self.work_dir.mkdir(exist_ok=True)

        # Get genome
        if search_term is not None:
            assert email is not None, "Email is required for NCBI API"

            # Hash the search term to use as filename in cache dir
            search_hashed = hashlib.sha1(search_term.encode()).hexdigest()
            genome_file = self.work_dir / f"db_{search_hashed}.gbk"
            genome_fasta = genome_file.with_suffix(".fasta")

            self._genome = download_genome(search_term, email, retmax, genome_file)
        else:
            genome_file = Path(genome_file)
            genome_fasta = self.work_dir / (genome_file.stem + ".fasta")

            if genome_file.suffix == ".gff":
                self._genome = SeqIO.to_dict(GFF.parse(genome_file))
            elif genome_file.suffix == ".gbk":
                self._genome = SeqIO.to_dict(SeqIO.parse(genome_file, "genbank"))
            else:
                raise RuntimeError(
                    "Wrong format expected either `.gbk` or `.gff` file."
                )

        # Save in fasta format (only acceptable format for makeblastdb)
        SeqIO.write(
            [x for x in self._genome.values() if x.seq.defined],
            genome_fasta,
            "fasta",
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

        self._seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

        self._seq_ids = tuple(
            {
                re.sub(f"{self._fwd_suf}$|{self._rev_suf}$", "", x)
                for x in self.seqs.keys()
            }
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

    def get_inserts(
        self,
        seq_id_or_idx,
        insert_types="both",
        filter_threshold=None,
        insert_max_len=10000,
    ):
        if isinstance(seq_id_or_idx, int):
            seq_id = self._seq_ids[seq_id_or_idx]
        else:
            seq_id = seq_id_or_idx

        matched = []

        fwds, revs = [], []
        if seq_id + self._fwd_suf in self._blast_results:
            fwds = self._blast_results[seq_id + self._fwd_suf].hsps
        if seq_id + self._rev_suf in self._blast_results:
            revs = self._blast_results[seq_id + self._rev_suf].hsps

        if filter_threshold is not None:
            if not 0 <= filter_threshold <= 1:
                raise ValueError("Filter value needs to be between 0 and 1")
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

        matched = [
            Insert(
                *x,
                seq_id=seq_id,
                seq_len=[
                    len(self._seqs[x[0].query_id]),
                    len(self._seqs[x[1].query_id]),
                ],
            )
            for x in matched
        ]
        unmatched = [
            Insert(
                x,
                seq_id=seq_id,
                seq_len=[
                    len(self._seqs[x.query_id]),
                ],
            )
            for x in unmatched
        ]

        if insert_types == "matched":
            return matched
        elif insert_types == "unmatched":
            return unmatched
        else:
            return matched + unmatched

    def get_all_inserts(
        self, insert_types="both", filter_threshold=None, insert_max_len=10000
    ):
        return [
            x
            for seq_id in self.seq_ids
            for x in self.get_inserts(
                seq_id,
                insert_types=insert_types,
                filter_threshold=filter_threshold,
                insert_max_len=insert_max_len,
            )
        ]

    def to_dataframe(self, insert_types="both", filter_threshold=None):
        return pd.DataFrame(
            [
                (seq_id, x.hit_id, x.start, x.end, x.strand, len(x))
                for seq_id in self.seq_ids
                for x in self.get_inserts(
                    seq_id, insert_types=insert_types, filter_threshold=filter_threshold
                )
            ],
            columns=(
                "seq_id",
                "NCBI_accession_number",
                "insert_start",
                "insert_end",
                "insert_strand",
                "insert_length",
            ),
        )

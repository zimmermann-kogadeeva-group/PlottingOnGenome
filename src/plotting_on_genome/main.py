#!/usr/bin/env python3

import hashlib
import json
import re
import subprocess
from copy import deepcopy
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


def get_endpoints(hsp1, hsp2=None, avg_insert_length=4000):
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


# TODO: use dataclass
class Insert(object):
    def __init__(self, start, end, seq_id, cov, hsp1, hsp2=None, genes=None):
        # TODO: maybe switch to using a list of hsps instead of hsp1 and hsp2
        self.hsp1 = hsp1
        self.hsp2 = hsp2
        self.strand = self.hsp1.hit_strand
        self.hit_id = self.hsp1.hit_id
        self.query_id = self.hsp1.query_id
        self.seq_id = seq_id
        self.coverage = cov
        self.start = start
        self.end = end
        self.genes = genes

    def __len__(self):
        return self.end - self.start

    def to_dataframe(self):
        return pd.DataFrame(
            [
                {
                    "start": gene.location.start,
                    "end": gene.location.end,
                    "strand": gene.location.strand,
                    "type": gene.type,
                    **gene.qualifiers,
                }
                for gene in self.genes
            ]
        )


class Pipeline(object):
    def __init__(
        self,
        seq_file,
        work_dir,
        *,
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
        self.work_dir.mkdir(exist_ok=True, parents=True)

        # Get genome
        if search_term is not None:
            assert email is not None, "Email is required for NCBI API"

            # Hash the search term to use as filename in cache dir
            search_hashed = hashlib.sha1(search_term.encode()).hexdigest()
            genome_file = self.work_dir / f"db_{search_hashed}.gbk"
            genome_fasta = genome_file.with_suffix(".fasta")

            # TODO: use diskcache package instead of hashing it yourself
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

    def get_genes(self, start, end, hit_id, buffer=4000):
        start_ = start - buffer
        end_ = end + buffer

        return [
            shift_feature(gene, start_)
            for gene in self._genome[hit_id][start_:end_].features
        ]

    def get_inserts(
        self,
        seq_id_or_idx,
        insert_types="both",
        filter_threshold=None,
        insert_max_len=10000,
        buffer=4000,
        avg_insert_length=4000,
    ):
        assert insert_types in ("matched", "unmatched", "both")
        if isinstance(seq_id_or_idx, int):
            seq_id = self._seq_ids[seq_id_or_idx]
        else:
            seq_id = seq_id_or_idx

        # Get all relevant fwd and rev hits from BLAST
        fwds, revs = [], []
        if seq_id + self._fwd_suf in self._blast_results:
            fwds = self._blast_results[seq_id + self._fwd_suf].hsps
        if seq_id + self._rev_suf in self._blast_results:
            revs = self._blast_results[seq_id + self._rev_suf].hsps

        # Apply the coverage filter
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

        # Match appropriate hits
        matched = []
        unmatched = set(fwds).union(revs)
        for fwd, rev in product(fwds, revs):
            if (
                fwd.hit_id == rev.hit_id
                and fwd.hit_strand == -rev.hit_strand
                and 0 < get_length(fwd, rev) < insert_max_len
            ):
                matched.append([fwd, rev])
                unmatched -= {fwd, rev}

        # TODO: improve below
        new_matched = []
        for fwd, rev in matched:
            start, end = get_endpoints(fwd, rev)
            genes = self.get_genes(start, end, fwd.hit_id, buffer)
            cov = [
                len(fwd.query.seq) / len(self.seqs[fwd.query_id]),
                len(rev.query.seq) / len(self.seqs[rev.query_id]),
            ]

            new_matched.append(Insert(start, end, seq_id, cov, fwd, rev, genes=genes))

        new_unmatched = []
        for seq in unmatched:
            # TODO: maybe move avg_insert_length to object attribute
            start, end = get_endpoints(seq, avg_insert_length=avg_insert_length)
            genes = self.get_genes(start, end, seq.hit_id, buffer)
            cov = [len(seq.query.seq) / len(self.seqs[seq.query_id])]

            new_unmatched.append(Insert(start, end, seq_id, cov, seq, genes=genes))

        if insert_types == "matched":
            return new_matched
        elif insert_types == "unmatched":
            return new_unmatched
        else:
            return new_matched + new_unmatched

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

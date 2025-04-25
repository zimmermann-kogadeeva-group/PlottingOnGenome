import hashlib
import subprocess
from copy import deepcopy
from pathlib import Path

import pandas as pd
from BCBio import GFF
from Bio import Entrez, SearchIO, SeqIO


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def _download_genome(search_term, email, retmax=100):
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

    return data_gb


def download_genome(search_term, email, retmax=100, output_path=None):

    if output_path is None:
        data_gb = _download_genome(search_term, email, retmax)

    else:
        output_path = Path(output_path)

        # Check output_path is given
        if Path(output_path).exists():
            # Read in the file
            data_gb = SeqIO.to_dict(SeqIO.parse(output_path, "genbank"))

        # Otherwise retrieve the records from NCBI
        else:
            data_gb = _download_genome(search_term, email, retmax)
            SeqIO.write(data_gb.values(), output_path, "genbank")

    return data_gb


def get_genome_file(work_dir, search_term, retmax):
    work_dir = Path(work_dir)
    # Hash the search term to use as filename in cache dir
    search_hashed = hashlib.sha1((search_term + str(retmax)).encode()).hexdigest()
    return work_dir / f"db_{search_hashed}.gbk"


def get_genome(genome_file, genome_fasta, search_term=None, email=None, retmax=None):
    # Get genome
    if search_term is not None:
        if email is None:
            raise ValueError("Email is required for NCBI API")

        genome = download_genome(search_term, email, retmax, genome_file)
    else:
        if genome_file.suffix == ".gff":
            genome = SeqIO.to_dict(GFF.parse(genome_file))
        elif genome_file.suffix == ".gbk":
            genome = SeqIO.to_dict(SeqIO.parse(genome_file, "genbank"))
        else:
            raise RuntimeError("Wrong format expected either `.gbk` or `.gff` file.")

    defined_seqs = [x for x in genome.values() if x.seq.defined]

    if not len(defined_seqs):

        message = f"Genome file {genome_file} does not contain any sequences!"
        if search_term is not None:
            message = (
                f"Search term '{search_term}' did not yield a genome with sequences"
            )
        raise RuntimeError(message)

    # Save in fasta format (only acceptable format for makeblastdb)
    SeqIO.write(defined_seqs, genome_fasta, "fasta")

    return genome


def _correct_hit_id(x):
    """Helper function to get rid red| or emb| added by BLAST to contig ID"""
    if "|" in x.id:
        x.id = x.id.split("|")[1]
    return x


def run_blast(seq_file, db_file, blast_output):
    # TODO: check how to catch errors from blast
    # make blast database
    run1 = subprocess.run(
        f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run1.returncode != 0:
        raise OSError(f"Failed to run makeblastdb: {run1.stderr.decode()}")

    # align input sequences with query using blastn with default options
    run2 = subprocess.run(
        f"blastn -query {seq_file} -db {db_file} -out {blast_output} -outfmt 5",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run2.returncode != 0:
        raise OSError(f"Failed to run blastn: {run2.stderr.decode()}")

    # Return a dictionary correcting the hit IDs which get an unnecessary
    # prefix from BLAST
    return {
        x.id: x.hit_map(_correct_hit_id)
        for x in SearchIO.parse(blast_output, "blast-xml")
    }


def get_inserts_df(inserts_dicts, insert_type="both", filter_threshold=None, **kwargs):

    inserts_dfs = [
        x.to_dataframe(
            insert_type=insert_type, filter_threshold=filter_threshold
        ).assign(genome=name)
        for name, x in inserts_dicts.items()
    ]

    return pd.concat(inserts_dfs, ignore_index=True)


def get_insert_presence_df(
    insert_dicts, insert_type="both", filter_threshold=None, **kwargs
):
    dfs = [
        pd.DataFrame(
            x.get_insert_ids(insert_type, filter_threshold),
            columns=["insert_ids"],
        ).assign(genome=name)
        for name, x in insert_dicts.items()
    ]

    df_insert_presence = (
        pd.concat(dfs, ignore_index=True)
        .groupby(["insert_ids", "genome"], as_index=False)
        .agg(num_inserts=pd.NamedAgg("insert_ids", "count"))
        .pivot(index="insert_ids", columns="genome")
        .loc[:, ("num_inserts")]
    )
    return df_insert_presence


def get_genes_df(
    inserts_dicts, insert_type="both", filter_threshold=None, buffer=None, **kwargs
):
    genes_dfs = [
        x.genes_to_dataframe(
            insert_type=insert_type,
            filter_threshold=filter_threshold,
            buffer=buffer,
        ).assign(genome=name)
        for name, x in inserts_dicts.items()
    ]

    return pd.concat(genes_dfs, ignore_index=True)

import subprocess
from copy import deepcopy
from pathlib import Path

from Bio import Entrez, SearchIO, SeqIO


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


def _correct_hit_id(x):
    """Helper function to get rid red| or emb| added by BLAST to contig ID"""
    if "|" in x.id:
        x.id = x.id.split("|")[1]
    return x


def run_blast(seq_file, db_file, blast_output):
    # TODO: check how to catch errors from blast
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
        x.id: x.hit_map(_correct_hit_id)
        for x in SearchIO.parse(blast_output, "blast-xml")
    }

from pathlib import Path
from tempfile import TemporaryDirectory

import streamlit as st

import plotting_on_genome as pog


class TempDirManager:
    def __init__(self, dirpath=None):
        self.dirpath = dirpath

    def __enter__(self):
        if self.dirpath is None:
            self.temp_dir = TemporaryDirectory()
            return self.temp_dir.name
        else:
            Path(self.dirpath).mkdir(exist_ok=True)
            return self.dirpath

    def __exit__(self, exc_type, exc_value, traceback):
        if self.dirpath is None:
            self.temp_dir.cleanup()


def get_main_inputs(workdir=False):

    # Radio buttons to switch between text input and file upload
    genome_src = st.radio("Genome:", ("NCBI", "File"))

    # Text inputs
    genome_fh = None
    search_terms = None
    email = None
    retmax = None

    if genome_src == "NCBI":

        search_terms = st.text_input(
            "Search term for genome (separete different search terms with a comma):",
            key="search_term",
        ).split(",")

        email = st.text_input("Email:", key="email")
        retmax = st.text_input(
            "Retmax (Max. # of records from NCBI)", 200, key="retmax"
        )
    else:
        genome_fh = st.file_uploader(
            "Upload genome:",
            type=("gbk", "gff"),
            key="genome",
            accept_multiple_files=True,
        )

    seq_fh = st.file_uploader(
        "Sequence file:",
        type=["fasta", "fas", "fna"],
        key="seqs",
        accept_multiple_files=True,
    )

    fwd_suf, rev_suf = None, None
    if st.toggle("Match forward and reverse sequences"):
        fwd_suf = st.text_input("Forward suffix:", "_F", key="fwd_suf")
        rev_suf = st.text_input("Reverse suffix:", "_R", key="rev_suf")

    workdir_path = None
    if workdir:
        workdir_path = st.text_input("workdir", "Output")

    most_inputs = {
        "seq_fh": seq_fh,
        "email": email,
        "retmax": retmax,
        "fwd_suf": fwd_suf,
        "rev_suf": rev_suf,
        "workdir": workdir_path,
    }

    st.write(search_terms)

    all_inputs = []
    if search_terms is not None:
        all_inputs = {
            x: most_inputs.copy() | {"search_term": x, "genome_fh": None}
            for x in search_terms
        }
    else:
        all_inputs = {
            x.name: most_inputs.copy() | {"genome_fh": x, "search_term": None}
            for x in genome_fh
        }

    return all_inputs


@st.cache_data
def run_pipeline(
    seq_fh,
    genome_fh,
    email,
    search_term,
    retmax,
    fwd_suf,
    rev_suf,
    workdir=None,
):
    with TempDirManager(workdir) as work_dir:
        dirpath = Path(work_dir)
        if genome_fh is not None:
            genome_path = str(dirpath / genome_fh.name)
            with open(genome_path, "wb") as fh:
                fh.write(genome_fh.getvalue())
        else:
            genome_path = None

        seq_path = str(dirpath / "seqs.fasta")
        with open(seq_path, "wb") as fh:
            for seq in seq_fh:
                # Make sure that there is new line between individual seqs
                fh.write(seq.getvalue())
                fh.write(b"\n")

        try:
            res = pog.Mapping(
                seq_file=seq_path,
                work_dir=dirpath,
                genome_file=genome_path,
                search_term=search_term,
                email=email,
                retmax=retmax,
                fwd_suffix=fwd_suf,
                rev_suffix=rev_suf,
            )
        except RuntimeError as e:
            st.error(e)
            return None
        else:
            return res

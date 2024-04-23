#!/usr/bin/env python3

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import streamlit as st

from plotting_on_genome import Pipeline, plot_insert, plot_insert_dists, plot_on_genome

# st.set_page_config(layout="wide")

if "stage" not in st.session_state:
    st.session_state.stage = 0

if "filter_threshold" not in st.session_state:
    st.session_state.filter_threshold = None

if "insert_types" not in st.session_state:
    st.session_state.insert_types = "both"

if "workdir" not in st.session_state:
    st.session_state.workdir = None


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


@st.cache_data
def run_pipeline(
    seq_fh, genome_fh, email, search_term, retmax, fwd_suf, rev_suf, workdir=None
):
    with st.spinner("Processing..."), TempDirManager(workdir) as work_dir:
        dirpath = Path(work_dir)
        if genome_fh is not None:
            genome_path = str(dirpath / genome_fh.name)
            with open(genome_path, "wb") as fh:
                fh.write(genome_fh.getbuffer())
        else:
            genome_path = None

        seq_path = str(dirpath / seq_fh.name)
        with open(seq_path, "wb") as fh:
            fh.write(seq_fh.getbuffer())

        return Pipeline(
            seq_path,
            dirpath,
            genome_path,
            email,
            search_term,
            retmax,
            fwd_suf,
            rev_suf,
        )


def submit(*args):
    st.session_state.pipeline = run_pipeline(*args)
    st.session_state.stage = 1


def get_main_inputs():
    parser = argparse.ArgumentParser("PlottingOnGenome")
    parser.add_argument("--workdir", action="store_true")
    args = parser.parse_args()

    # Radio buttons to switch between text input and file upload
    genome_src = st.radio("Genome:", ("NCBI", "File"))

    # Text inputs
    genome_fh = None
    search_term = None
    email = None
    retmax = None

    if genome_src == "NCBI":
        search_term = st.text_input("Search term:", key="search_term")
        email = st.text_input("Email:", key="email")
        retmax = st.text_input(
            "Retmax (Max. # of records from NCBI)", 200, key="retmax"
        )
    else:
        genome_fh = st.file_uploader("Upload genome:", type="gbk", key="genome")

    seq_fh = st.file_uploader("Sequence file:", type=["fasta", "fna"], key="seqs")

    fwd_suf = st.text_input("Forward suffix:", "_F", key="fwd_suf")
    rev_suf = st.text_input("Reverse suffix:", "_R", key="rev_suf")

    filter_threshold = st.text_input("Filter threshold:", None)
    if filter_threshold is not None:
        st.session_state.filter_threshold = int(filter_threshold)

    st.session_state.insert_types = st.selectbox(
        "insert types:", ["both", "matched", "unmatched"]
    )

    if args.workdir:
        st.session_state.workdir = st.text_input("workdir", "Output")

    st.button(
        "Submit",
        on_click=submit,
        args=[
            seq_fh,
            genome_fh,
            email,
            search_term,
            retmax,
            fwd_suf,
            rev_suf,
            st.session_state.workdir,
        ],
    )


def main():
    st.title("PlottingOnGenome")

    if st.session_state.stage == 0:
        get_main_inputs()

    if st.session_state.stage >= 1:
        option = st.selectbox(
            "plot type:",
            ["plot inserts", "plot all on genome", "plot insert dists", "show data"],
            None,
        )
        p = st.session_state.pipeline
        insert_types = st.session_state.insert_types
        filter_threshold = st.session_state.filter_threshold

        if option == "plot inserts":
            seq_id = st.selectbox("Select sequence id:", p.seq_ids)

            inserts = p.get_inserts(seq_id, insert_types, filter_threshold)

            for insert in inserts:
                fig, axs = plt.subplots(2, 1, figsize=(10, 10), height_ratios=[3, 5])
                axs = plot_insert(insert, p.genome, axs=axs)
                st.pyplot(fig)
                plt.close()

        elif option == "plot all on genome":
            labels = st.radio("Labels:", (True, False))
            inserts = p.get_all_inserts(insert_types, filter_threshold)

            fig, ax = plt.subplots(figsize=(10, 10 * (1 + 2 * labels)))
            ax = plot_on_genome(inserts, p.genome, labels=labels, ax=ax)
            st.pyplot(fig, use_container_width=True)

        elif option == "plot insert dists":
            inserts = p.get_all_inserts(insert_types, filter_threshold)

            fig, axs = plt.subplots(1, 3, figsize=(12, 5))
            plot_insert_dists(inserts, axs)
            st.pyplot(fig)

        elif option == "show data":
            st.write(p.to_dataframe(insert_types, filter_threshold))

        if st.button("Reset"):
            st.session_state.stage = 0
            st.rerun()


if __name__ == "__main__":
    main()

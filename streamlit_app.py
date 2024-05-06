#!/usr/bin/env python3

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import plotting_on_genome as pog

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

        seq_path = str(dirpath / "seqs.fasta")
        with open(seq_path, "wb") as fh:
            for seq in seq_fh:
                fh.write(seq.getbuffer())

        res = pog.Pipeline(
            seq_file=seq_path,
            work_dir=dirpath,
            genome_file=genome_path,
            search_term=search_term,
            email=email,
            retmax=retmax,
            fwd_suffix=fwd_suf,
            rev_suffix=rev_suf,
        )

        return res


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
        genome_fh = st.file_uploader(
            "Upload genome:", type=("gbk", "gff"), key="genome"
        )

    seq_fh = st.file_uploader(
        "Sequence file:",
        type=["fasta", "fas", "fna"],
        key="seqs",
        accept_multiple_files=True,
    )

    fwd_suf = st.text_input("Forward suffix:", "_F", key="fwd_suf")
    rev_suf = st.text_input("Reverse suffix:", "_R", key="rev_suf")

    filter_threshold = st.text_input("Filter threshold (optional):", None)
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


def show_results():
    option = st.selectbox(
        "plot type:",
        [
            "plot inserts",
            "plot genome",
            "plot insert dists",
        ],
        None,
    )
    p = st.session_state.pipeline
    insert_types = st.session_state.insert_types
    filter_threshold = st.session_state.filter_threshold

    if p is not None:
        if option == "plot inserts":
            seq_id = st.selectbox("Select sequence id:", p.seq_ids)

            inserts = p.get_inserts(seq_id, insert_types, filter_threshold)

            st.write(
                pd.concat(
                    [
                        insert.to_dataframe().assign(insert_idx=idx + 1)
                        for idx, insert in enumerate(inserts)
                    ]
                )
            )
            for idx, insert in enumerate(inserts):
                fig, axs = plt.subplots(2, 1, figsize=(10, 10), height_ratios=[3, 5])
                fig.suptitle(f"Insert {idx+1}")
                axs = pog.plot_insert(insert, axs=axs)
                st.pyplot(fig)
                plt.close()

        elif option == "plot genome":
            labels = st.toggle("Labels")
            inserts = []
            if st.toggle("Plot inserts", True):
                inserts = p.get_all_inserts(insert_types, filter_threshold)

            st.write(p.to_dataframe(insert_types, filter_threshold))

            fig, ax = plt.subplots(figsize=(10, 10 * (1 + 2 * labels)))
            ax = pog.plot_on_genome(p.genome, inserts, labels=labels, ax=ax)
            st.pyplot(fig, use_container_width=True)
            plt.close()

        elif option == "plot insert dists":
            plot_type = st.radio(
                "Plot type:", ["histogram", "violinplot+boxplot+stripplot"]
            )
            inserts = p.get_all_inserts(insert_types, filter_threshold)

            if plot_type == "histogram":
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))
                pog.plot_histogram(inserts, axs)
                st.pyplot(fig)
            elif plot_type == "violinplot+boxplot+stripplot":
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))
                pog.plot_dists(inserts, axs)
                st.pyplot(fig)
            else:
                raise ValueError("Incorrect plot type")

            plt.close()

        if st.button("Reset"):
            st.session_state.stage = 0
            st.rerun()


def main():
    st.title("PlottingOnGenome")

    if st.session_state.stage == 0:
        get_main_inputs()

    if st.session_state.stage >= 1:
        show_results()


if __name__ == "__main__":
    main()

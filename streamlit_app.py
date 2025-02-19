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
                fh.write(genome_fh.getvalue())
        else:
            genome_path = None

        seq_path = str(dirpath / "seqs.fasta")
        with open(seq_path, "wb") as fh:
            for seq in seq_fh:
                # Make sure that there is new line between individual seqs
                fh.write(seq.getvalue())
                fh.write(b"\n")

        res = pog.InsertsDict(
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
        ["plot inserts", "plot genome", "plot insert dists"],
        None,
    )
    inserts_all = st.session_state.pipeline
    insert_types = st.session_state.insert_types

    filter_threshold = st.slider("Filter threshold:", 0.0, 1.0, 0.7)

    if inserts_all is not None:
        if option == "plot inserts":

            buffer = st.slider(
                "View window size:",
                min_value=0,
                max_value=10000,
                value=4000,
                help="Number of bases either side of the insert",
            )

            st.write("Features to dispaly:")
            ft_checkboxes = {
                "CDS": st.checkbox("CDS", value=True),
                "gene": st.checkbox("gene", value=True),
            }
            feature_types = [k for k, v in ft_checkboxes.items() if v]

            seq_id = st.selectbox("Select sequence id:", inserts_all.seq_ids)

            inserts = inserts_all.get(seq_id, insert_types, filter_threshold)

            if len(inserts):

                df_features = pd.concat(
                    [
                        insert.to_dataframe().assign(insert_idx=idx + 1)
                        for idx, insert in enumerate(inserts)
                    ]
                ).reset_index(drop=True)

                st.write(df_features)

                for idx, insert in enumerate(inserts):
                    fig, axs = plt.subplots(
                        2, 1, figsize=(10, 10), height_ratios=[3, 5]
                    )
                    fig.suptitle(f"Insert {idx+1}")
                    axs = insert.plot(
                        buffer=buffer, axs=axs, feature_types=feature_types
                    )
                    st.pyplot(fig)
                    plt.close()
            else:
                st.write(f"No inserts found for {seq_id}!")

        elif option == "plot genome":
            labels = st.toggle("Labels")
            inserts = []
            if st.toggle("Plot inserts", True):
                inserts = inserts_all.filter(insert_types, filter_threshold)

            st.write(inserts_all.to_dataframe(insert_types, filter_threshold))

            fig, ax = plt.subplots(figsize=(10, 10 * (1 + 2 * labels)))
            ax = inserts_all.plot(inserts, show_labels=labels, ax=ax)
            st.pyplot(fig, use_container_width=True)
            plt.close()

        elif option == "plot insert dists":
            plot_type = st.radio(
                "Plot type:", ["histogram", "violinplot+boxplot+stripplot"]
            )
            inserts = inserts_all.filter(insert_types, filter_threshold)

            if plot_type == "histogram":
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))
                pog.plot_histogram(inserts, axs=axs)
                st.pyplot(fig)

            elif plot_type == "violinplot+boxplot+stripplot":
                fig, axs = plt.subplots(1, 2, figsize=(12, 5))
                pog.plot_dists(inserts, axs=axs)
                st.pyplot(fig)

            else:
                raise ValueError("Incorrect plot type")

            plt.close()


def main():
    col1, col2 = st.columns((6, 1), gap="large", vertical_alignment="bottom")
    with col1:
        st.title("PlottingOnGenome")
    with col2:
        if st.session_state.stage == 1:
            if st.button("Reset"):
                st.session_state.stage = 0
                st.rerun()

    if st.session_state.stage == 0:
        get_main_inputs()

    if st.session_state.stage >= 1:
        show_results()


if __name__ == "__main__":
    main()

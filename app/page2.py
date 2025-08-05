from itertools import chain

import pandas as pd
import streamlit as st
from plotting import plot_genomes, plot_inserts, plot_inserts_dist

import plotting_on_genome as pog


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def sidebar_opts():
    if st.sidebar.button("Reset", use_container_width=True):
        st.session_state.stage = 0
        st.session_state.search_term_count = 1
        st.rerun()

    insert_type = st.sidebar.selectbox(
        "Insert type:",
        ["both", "paired", "unpaired"],
        help=(
            "There are two types of inserts: paired and unpaired. Paired are when "
            "both forward and reverse sequences could be paired to one another. "
            "Unpaired are cases when either forward or reverse sequence could not be "
            "paired to the corresponding one."
        ),
    )

    filter_threshold = st.sidebar.slider(
        "Filter threshold:",
        min_value=0.0,
        max_value=1.0,
        value=0.7,
        help=(
            "Threshold value for insert coverage value. Insert coverage here is the "
            "ratio of length of sequence mapped to genome to the length of the whole "
            "sequence."
        ),
    )

    buffer = st.sidebar.slider(
        "View window size:",
        min_value=0,
        max_value=10000,
        value=4000,
        help="Number of bases either side of the insert",
    )

    params = dict(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )

    return params


def download_tables(df_inserts, df_genes):

    st.write(
        df_inserts.groupby(["genome", "insert_paired"], as_index=False)
        .agg(num_inserts=pd.NamedAgg("seq_id", "nunique"))
        .pivot(index="insert_paired", columns="genome", values="num_inserts")
        .reindex(index=[False, True])
        .rename(index={False: "Unpaired", True: "Paired"})
        .rename_axis(index="")
        .fillna(0)
    )

    st.download_button(
        label="Download all mapped seqs as CSV",
        data=df_inserts.pipe(convert_df, index=False),
        file_name="inserts.csv",
        mime="text/csv",
        use_container_width=True,
    )

    st.download_button(
        label="Download all genes as CSV",
        data=df_genes.pipe(convert_df, index=False),
        file_name="genes.csv",
        mime="text/csv",
        use_container_width=True,
    )


def select_genomes(genome_labels):
    st.sidebar.write("Genomes:")
    res_choice = [
        name
        for idx, name in enumerate(genome_labels)
        if st.sidebar.checkbox(name, True)
    ]

    return res_choice


def select_seq_id(possible_seq_ids, clusters_dict, st_key):

    sel = st.multiselect(
        "Select sequence id:",
        possible_seq_ids + list(clusters_dict),
        None,
        format_func=lambda x: (
            f"{x[0]} - cluster {x[1]}" if isinstance(x, (list, tuple)) else x
        ),
        key=st_key,
    )

    sel = set(sel)
    clusters = {x: clusters_dict[x] for x in sel if x in clusters_dict}
    seq_id = list(sel - set(clusters.keys()))

    return seq_id, clusters


def get_inserts_cond(seq_ids, clusters):
    cluster_ins_ids = list(chain.from_iterable(clusters.values()))
    return " or ".join(
        [
            (
                f"(seq_id == '{seq_id[0]}' and  insert_idx == {seq_id[1]})"
                if len(seq_id) == 2
                else f"seq_id == '{seq_id}'"
            )
            for seq_id in set(seq_ids + cluster_ins_ids)
        ]
    )


def show_results():

    if st.session_state.results is not None:
        all_results = pog.Comparison(st.session_state.results)

        params = sidebar_opts()

        tabbed = st.sidebar.toggle("Tabbed interface", True)
        if tabbed:
            genome_view, insert_view = st.tabs(["Genome view", "Insert view"])
        else:
            genome_view, insert_view = st.columns(2)

        res_choice = select_genomes(all_results.keys())

        if len(res_choice):
            df_insert_presence = all_results.get_insert_presence_df(
                res_choice, **params
            )
            df_inserts = all_results.get_inserts_df(res_choice, **params)
            df_genes = all_results.get_genes_df(res_choice, **params).map(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )

            # Possible seq_ids and clusters
            pos_seq_ids = df_insert_presence.index.tolist()
            pos_clusters = all_results.get_clusters(res_choice, **params)

            seq_id = []
            with genome_view:
                st.header("Genome view")

                # download all genes and inserts
                download_tables(df_inserts, df_genes)

                with st.expander("Mapped seqs presence/absence table"):
                    st.write(df_insert_presence)

                # Select inserts
                seq_id, clust_sel = select_seq_id(pos_seq_ids, pos_clusters, "seq1")

                if seq_id or clust_sel:
                    df_inserts = df_inserts.query(get_inserts_cond(seq_id, clust_sel))

                with st.expander("Table of mapped seqs"):
                    st.write(df_inserts)

                with st.expander("Distribution of mapped seqs"):
                    plot_inserts_dist(df_inserts)

                plot_genomes(
                    all_results,
                    res_choice,
                    seq_id or pos_seq_ids,
                    **params,
                )

            with insert_view:
                st.header("Insert view")
                if tabbed:
                    seq_id, clust_sel = select_seq_id(pos_seq_ids, pos_clusters, "seq2")

                genome_choice = st.multiselect("Genome:", res_choice, None)
                if not len(genome_choice):
                    genome_choice = res_choice

                plot_inserts(all_results, genome_choice, seq_id, clust_sel, **params)

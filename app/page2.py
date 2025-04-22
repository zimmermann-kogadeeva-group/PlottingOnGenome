from collections import defaultdict

import matplotlib.pyplot as plt
import streamlit as st

import plotting_on_genome as pog


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def download_tables(all_inserts, insert_type, filter_threshold, buffer, **kwargs):
    df_inserts = all_inserts.to_dataframe(
        insert_type=insert_type, filter_threshold=filter_threshold
    )
    df_genes = all_inserts.genes_to_dataframe(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )

    counts = defaultdict(
        int, df_inserts.groupby("insert_matched").insert_matched.count().to_dict()
    )
    st.write(
        f"Number of inserts:\n "
        f"- matched: {counts[True]}\n "
        f"- unmatched: {counts[False]}"
    )

    st.download_button(
        label="Download all inserts as CSV",
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


def plot_inserts(all_inserts, insert_type, filter_threshold, buffer, **kwargs):
    st.write("Features to display:")
    ft_checkboxes = {
        "CDS": st.checkbox("CDS", value=True),
        "gene": st.checkbox("gene", value=True),
    }
    feature_types = [k for k, v in ft_checkboxes.items() if v]

    colorbar = st.toggle("Color genes by overlap")

    seq_id = st.selectbox("Select sequence id:", all_inserts.seq_ids)

    inserts = all_inserts.get(
        seq_id, insert_type=insert_type, filter_threshold=filter_threshold
    )

    df_genes = all_inserts.genes_to_dataframe(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )

    if len(inserts):

        st.write(df_genes.query("seq_id == @seq_id"))

        for insert in inserts:
            st.write(
                f"Insert {insert.idx}: "
                f"coverage = {insert.coverage:.2f}, "
                f"matched = {insert.matched}"
            )
            fig, axs = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[3, 5])
            fig.suptitle(f"Insert {insert.idx}")
            axs = insert.plot(
                buffer=buffer,
                axs=axs,
                feature_types=feature_types,
                colorbar=colorbar,
            )
            st.pyplot(fig)
            plt.close()
    else:
        st.write(f"No inserts found for {seq_id}!")


def plot_genome(all_inserts, insert_type, filter_threshold, **kwargs):

    labels = st.toggle("Labels")
    insert_idxs = None if st.toggle("Plot inserts", True) else []

    df_inserts = all_inserts.to_dataframe(
        insert_type=insert_type, filter_threshold=filter_threshold
    )
    st.write(df_inserts)

    fig, ax = plt.subplots(figsize=(10, 10 * (1 + 2 * labels)))
    ax = all_inserts.plot(
        insert_idxs,
        show_labels=labels,
        ax=ax,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
    )
    st.pyplot(fig, use_container_width=True)
    plt.close()


def plot_dists(all_inserts, insert_type, filter_threshold, **kwargs):
    plot_type = st.radio("Plot type:", ["histogram", "violinplot+boxplot+stripplot"])
    inserts = all_inserts.get(
        insert_type=insert_type, filter_threshold=filter_threshold
    )

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


def show_results():
    inserts_all = st.session_state.pipeline

    insert_type = st.selectbox(
        "Insert type:",
        ["both", "matched", "unmatched"],
        help=(
            "There are two types of inserts: matched and unmatched. Matched are when "
            "both forward and reverse sequences could be matched to one another. "
            "Unmatched are cases when either forward or reverse sequence could not be "
            "matched to the corresponding one."
        ),
    )

    filter_threshold = st.slider(
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

    buffer = st.slider(
        "View window size:",
        min_value=0,
        max_value=10000,
        value=4000,
        help="Number of bases either side of the insert",
    )
    params = dict(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )

    if inserts_all is not None:
        # download all genes and inserts
        download_tables(inserts_all, **params)

        option = st.selectbox(
            "plot type:",
            ["plot inserts", "plot genome", "plot insert dists"],
            None,
        )

        if option == "plot inserts":
            plot_inserts(inserts_all, **params)

        elif option == "plot genome":
            plot_genome(inserts_all, **params)

        elif option == "plot insert dists":
            plot_dists(inserts_all, **params)

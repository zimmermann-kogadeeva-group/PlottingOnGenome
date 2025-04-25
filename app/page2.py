import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import plotting_on_genome as pog


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def download_tables(all_results, insert_type, filter_threshold, buffer, **kwargs):
    df_inserts = pog.get_inserts_df(all_results, insert_type, filter_threshold)

    df_genes = pog.get_genes_df(all_results, insert_type, filter_threshold, buffer).map(
        lambda x: ",".join(x) if isinstance(x, list) else x
    )

    st.write(
        df_inserts.groupby(["genome", "insert_matched"], as_index=False)
        .agg(num_inserts=pd.NamedAgg("insert_matched", "count"))
        .pivot(index="insert_matched", columns="genome", values="num_inserts")
        .reindex(index=[False, True])
        .rename(index={False: "Unmatched", True: "Matched"})
        .rename_axis(index="")
        .fillna(0)
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

    seq_id = st.selectbox(
        "Select sequence id:", all_inserts.get_insert_ids(insert_type, filter_threshold)
    )

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
    insert_type = st.sidebar.selectbox(
        "Insert type:",
        ["both", "matched", "unmatched"],
        help=(
            "There are two types of inserts: matched and unmatched. Matched are when "
            "both forward and reverse sequences could be matched to one another. "
            "Unmatched are cases when either forward or reverse sequence could not be "
            "matched to the corresponding one."
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

    all_results = st.session_state.results

    if all_results is not None:

        res_choice = [
            name
            for idx, name in enumerate(all_results.keys())
            if st.sidebar.checkbox(name, idx == 0)
        ]

        # download all genes and inserts
        download_tables({name: all_results[name] for name in res_choice}, **params)

        if len(res_choice) == 1:
            inserts_all = all_results[res_choice[0]]

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

        else:
            # TODO: add choice of seq_id
            res_subset = {name: all_results[name] for name in res_choice}

            df_insert_presence = pog.get_insert_presence_df(res_subset, **params)
            with st.expander("Insert presence/absence"):
                st.write(df_insert_presence)

            seq_id = st.selectbox("Select sequence id:", df_insert_presence.index, None)

            df_inserts = pog.get_inserts_df(res_subset, **params)
            if seq_id is not None:
                df_inserts = df_inserts.query("seq_id == @seq_id")

            with st.expander("Inserts info"):
                st.write(df_inserts)

            fig, ax = plt.subplots(figsize=(10, 10))
            ax = pog.plot_multiple_genomes(*res_subset.values(), seq_id=seq_id, ax=ax)
            st.pyplot(fig, use_container_width=True)
            plt.close()

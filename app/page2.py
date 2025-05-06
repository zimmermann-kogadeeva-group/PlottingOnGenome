import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st

import plotting_on_genome as pog


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def download_tables(df_inserts, df_genes):

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


def plot_inserts(all_inserts, seq_id, insert_type, filter_threshold, buffer, **kwargs):

    inserts = all_inserts.get(
        seq_id, insert_type=insert_type, filter_threshold=filter_threshold
    )

    # Get a table of genes and display it in the webapp
    df_genes = all_inserts.genes_to_dataframe(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )
    with st.expander("Genes table"):
        st.write(df_genes.query("seq_id == @seq_id"))

    feature_types = set()
    colorbar = False

    # Get display options from the user
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.checkbox("Display CDS", value=True):
            feature_types.update(("CDS",))
    with col2:
        if st.checkbox("Display genes", value=True):
            feature_types.update(("gene",))
    with col3:
        colorbar = st.checkbox("Color genes by overlap", value=False)

    if len(inserts):
        for insert in inserts:
            with st.expander(
                f"Insert {insert.idx}: "
                f"coverage = {insert.coverage:.2f}, "
                f"matched = {insert.matched}, "
                f"hit_id = {insert.hit_id}"
            ):
                fig, axs = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[2, 5])
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


def sidebar_opts():
    if st.sidebar.button("Reset", use_container_width=True):
        st.session_state.stage = 0
        st.session_state.search_term_count = 1
        st.rerun()

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

    return dict(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )


def plot_inserts_dist(data, palette="tab10"):
    if len(data) > 1:
        col1, col2 = st.columns(2)
        plot_type = "boxplot"
        sel_col = "insert_length"
        with col1:
            plot_type = st.radio("Plot type", ["boxplot", "violinplot", "scatterplot"])
        with col2:
            sel_col = st.selectbox("column", ["insert_length", "insert_coverage"])

        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.3)
        params = dict(
            data=data,
            x="genome",
            y=sel_col,
            hue="genome",
            ax=ax,
            palette=palette,
            hue_order=data.get("genome").unique(),
        )

        if plot_type == "violinplot":
            sns.violinplot(**params)
        if plot_type == "boxplot":
            sns.boxplot(**params)
        if plot_type == "scatterplot":
            sns.stripplot(**params, edgecolor="black", linewidth=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set(
            ylim=(
                0.5 * data.get(sel_col).min(),
                1.2 * data.get(sel_col).max(),
            )
        )

        st.pyplot(fig)
        plt.close()


def show_results():
    params = sidebar_opts()

    if st.session_state.results is not None:
        all_results = pog.Comparison(st.session_state.results)

        res_choice = [
            name
            for idx, name in enumerate(all_results.keys())
            if st.sidebar.checkbox(name, True)
        ]

        if len(res_choice):
            seq_id = None

            df_insert_presence = all_results.get_insert_presence_df(
                res_choice, **params
            )
            df_inserts = all_results.get_inserts_df(res_choice, **params)
            df_genes = all_results.get_genes_df(res_choice, **params).map(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )

            col1, col2 = st.columns(2)

            with col1:
                st.header("Genome view")
                # download all genes and inserts
                download_tables(df_inserts, df_genes)

                with st.expander("Insert presence/absence table"):
                    st.write(df_insert_presence)

                seq_id = st.selectbox(
                    "Select sequence id:", df_insert_presence.index, None
                )

                if seq_id is not None:
                    df_inserts = df_inserts.query("seq_id == @seq_id")

                with st.expander("Inserts info"):
                    st.write(df_inserts)

                with st.expander("Plot inserts dist"):
                    plot_inserts_dist(df_inserts)

                fig, ax = plt.subplots(figsize=(10, 10))
                ax = all_results.plot(res_choice, seq_id=seq_id, ax=ax, **params)
                st.pyplot(fig, use_container_width=True)
                plt.close()

            with col2:
                st.header("Insert view")
                if seq_id is not None:
                    genomes_list = (
                        df_insert_presence.loc[seq_id].dropna().index.tolist()
                    )

                    genome_choice = st.selectbox("Genome:", genomes_list, None)
                    if genome_choice is not None:
                        plot_inserts(all_results[genome_choice], seq_id, **params)

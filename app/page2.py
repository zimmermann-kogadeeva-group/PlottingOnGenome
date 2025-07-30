import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st

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

    return dict(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )


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
        plt.close(fig)


def plot_inserts(
    all_inserts, genome_choice, seq_id, insert_type, filter_threshold, buffer, **kwargs
):

    # Get a table of genes and display it in the webapp
    df_genes = all_inserts.get_genes_df(
        genome_choice,
        seq_id,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        buffer=buffer,
    )
    with st.expander("Genes table"):
        st.write(df_genes)

    feature_types = set()
    colorbar = False

    # Get display options from the user
    col1, col2 = "#ebf3ed", "#2e8b57"
    with st.expander("Plotting options"):
        if st.toggle("Display CDS", value=True):
            feature_types.update(("CDS",))
        if st.toggle("Display genes", value=True):
            feature_types.update(("gene",))
        colorbar = st.toggle("Color genes by overlap", value=False)
        if colorbar:
            col2 = st.color_picker("Inserts color:", value=col2)
            col1 = st.color_picker("Genes/CDS color:", value=col1)

    for genome in genome_choice:
        st.subheader(genome, divider=True)
        inserts = all_inserts[genome].get(
            seq_id,
            insert_type=insert_type,
            filter_threshold=filter_threshold,
        )
        for insert in inserts:
            st.write(
                f"Insert {insert.seq_id}: "
                f"index = {insert.idx}, "
                f"coverage = {insert.coverage:.2f}, "
                f"paired = {insert.paired}, "
                f"hit_id = {insert.hit_id}"
            )
            fig, axs = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[2, 5])
            fig.suptitle(f"Insert {insert.idx}")
            axs = insert.plot(
                buffer=buffer,
                axs=axs,
                feature_types=feature_types,
                colorbar=colorbar,
                col1=col1,
                col2=col2,
            )
            st.pyplot(fig)
            plt.close()


def plot_multiple_inserts(
    all_inserts, genome_choice, seq_id, insert_type, filter_threshold, buffer, **kwargs
):
    if isinstance(genome_choice, (tuple, list)):
        assert len(genome_choice) == 1, "Multiple genomes specified"
        genome_choice = genome_choice[0]

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

    fig, axs = plt.subplots(nrows=2, figsize=(10, 10))
    res = all_inserts[genome_choice]
    inserts = res.get(seq_id)
    for insert in inserts:
        st.write(
            f"Insert {insert.seq_id} with index {insert.idx}: "
            f"coverage = {insert.coverage:.2f}, "
            f"paired = {insert.paired}, "
            f"hit_id = {insert.hit_id}"
        )
    axs = res.plot_inserts(
        seq_id,
        axs=axs,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        buffer=buffer,
        feature_types=feature_types,
        colorbar=colorbar,
    )
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)


def plot_genomes(
    all_inserts, genome_choice, seq_id, insert_type, filter_threshold, buffer, **kwargs
):
    show_labels, show_titles = True, True
    num_cols = None

    with st.expander("Plotting options"):
        show_labels = st.toggle("Show contig labels", value=True)
        show_titles = st.toggle("Show genome labels", value=True)
        facet = st.toggle("Separate plot for each genome", value=False)
        if facet:
            num_cols = st.number_input("Number of columns", value=3)

    fig = all_inserts.plot(
        selection={g: seq_id for g in genome_choice},
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        show_labels=show_labels,
        show_titles=show_titles,
        facet_wrap=num_cols,
    )
    st.pyplot(fig, use_container_width=True)
    plt.close()


def get_inserts_cond(seq_ids):
    return " or ".join(
        [
            (
                f"(seq_id == '{seq_id[0]}' and  insert_idx == {seq_id[1]})"
                if len(seq_id) == 2
                else f"seq_id == '{seq_id}'"
            )
            for seq_id in seq_ids
        ]
    )


def show_results():
    params = sidebar_opts()

    tabbed = st.sidebar.toggle("Tabbed interface", True)
    if tabbed:
        genome_view, insert_view = st.tabs(["Genome view", "Insert view"])
    else:
        genome_view, insert_view = st.columns(2)

    if st.session_state.results is not None:
        all_results = pog.Comparison(st.session_state.results)

        st.sidebar.write("Genomes:")
        res_choice = [
            name
            for idx, name in enumerate(all_results.keys())
            if st.sidebar.checkbox(name, True)
        ]

        if len(res_choice):
            df_insert_presence = all_results.get_insert_presence_df(
                res_choice, **params
            )
            df_inserts = all_results.get_inserts_df(res_choice, **params)
            df_genes = all_results.get_genes_df(res_choice, **params).map(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )

            all_clusters = all_results.get_clusters(res_choice, **params)

            sel_seq_id = False
            seq_id = []
            with genome_view:
                st.header("Genome view")

                # download all genes and inserts
                download_tables(df_inserts, df_genes)

                with st.expander("Insert presence/absence table"):
                    st.write(df_insert_presence)

                # Select inserts
                seq_id = st.multiselect(
                    "Select sequence id:",
                    # [*insert_presence_df.index, *cluster_ids],
                    df_insert_presence.index,
                    None,
                )
                if not len(seq_id):
                    seq_id = df_insert_presence.index.tolist()
                    sel_seq_id = False
                else:
                    sel_seq_id = True

                if isinstance(seq_id, (tuple, list)) and len(seq_id):
                    df_inserts = df_inserts.query(get_inserts_cond(seq_id))

                with st.expander("Inserts info"):
                    st.write(df_inserts)

                with st.expander("Plot inserts dist"):
                    plot_inserts_dist(df_inserts)

                plot_genomes(all_results, res_choice, seq_id, **params)

            with insert_view:
                st.header("Insert view")
                if tabbed:
                    seq_id = st.multiselect(
                        "Select sequence id:",
                        # [*insert_presence_df.index, *cluster_ids],
                        df_insert_presence.index,
                        None,
                        key="seq_id_insert_view",
                    )
                    if not len(seq_id):
                        seq_id = df_insert_presence.index.tolist()
                        sel_seq_id = False
                    else:
                        sel_seq_id = True

                genomes_list = (
                    df_insert_presence.loc[seq_id, :]
                    .dropna(axis=1, how="all")
                    .columns.tolist()
                )

                if sel_seq_id:
                    genome_choice = st.multiselect("Genome:", genomes_list, None)
                    if not len(genome_choice):
                        genome_choice = res_choice

                    plot_inserts(all_results, genome_choice, seq_id, **params)

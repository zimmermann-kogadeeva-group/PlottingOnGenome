import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st


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
    all_inserts,
    genome_choice,
    seq_id,
    clusters=None,
    insert_type="both",
    filter_threshold=None,
    buffer=4000,
    **kwargs,
):

    # Get a table of genes and display it in the webapp
    df_genes = all_inserts.get_genes_df(
        genome_choice,
        seq_id,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        buffer=buffer,
    )
    with st.expander("Table of genes"):
        st.write(df_genes)

    # Set defaults for user-defined params
    feature_types = set()
    colorbar = False
    col1, col2 = "#ebf3ed", "#2e8b57"

    # Get display options from the user
    with st.expander("Plotting options"):
        if st.toggle("Display CDS", value=True):
            feature_types.update(("CDS",))
        if st.toggle("Display genes", value=True):
            feature_types.update(("gene",))

        colorbar = st.toggle("Color genes by overlap", value=False)

        col1 = st.color_picker("Genes/CDS color:", value=col1)
        col2 = st.color_picker("Inserts color:", value=col2)

    for genome in genome_choice:
        st.subheader(genome, divider="gray")
        inserts = all_inserts[genome].get(
            seq_id,
            insert_type=insert_type,
            filter_threshold=filter_threshold,
        )
        for insert in inserts:
            st.write(f"{insert:short}")
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

    if clusters is not None:
        for (genome, clust_idx), seq_ids in clusters.items():
            fig, axs = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[2, 5])
            fig.suptitle(f"Cluster {clust_idx}")
            clust_ins = all_inserts[genome].get(seq_ids, insert_type, filter_threshold)
            st.write(f"Cluster {clust_idx}:")
            st.write("\n".join([f"- {x:short}" for x in clust_ins]))
            axs = all_inserts[genome].plot_inserts(
                seq_ids,
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
    all_inserts, genome_choice, seq_id, insert_type, filter_threshold, **kwargs
):
    show_labels, show_titles = True, True
    num_cols = None

    with st.expander("Plotting options"):
        show_labels = st.toggle("Show contig labels", value=True)
        show_titles = st.toggle("Show genome labels", value=True)
        facet = st.toggle("Separate plot for each genome", value=False)
        if facet:
            num_cols = st.number_input("Number of columns", value=3, min_value=1)

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

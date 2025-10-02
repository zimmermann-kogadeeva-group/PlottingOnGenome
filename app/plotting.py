import zipfile
from collections import defaultdict
from io import BytesIO

import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from matplotlib import color_sequences


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
    comparison,
    genome_choices,
    seq_ids,
    clusters=None,
    insert_type="both",
    filter_threshold=None,
    buffer=4000,
    **kwargs,
):
    # Set defaults for user-defined params
    feature_types = set()
    colorbar = False
    col1, col2 = "#ebf3ed", "#2e8b57"

    clusters2 = clusters.other_repr

    # Get display options from the user
    with st.expander("Plotting options"):
        if st.toggle("Display CDS", value=True):
            feature_types.update(("CDS",))
        if st.toggle("Display genes", value=True):
            feature_types.update(("gene",))

        show_mappings = st.toggle(
            "Draw mapped sequence bounds",
            value=False,
            help=(
                "Draws additional vertical lines to indicate "
                "where fwd/rev sequences are mapped"
            ),
        )
        colorbar = st.toggle("Color genes by overlap", value=False)

        col1 = st.color_picker("Genes/CDS color:", value=col1)
        col2 = st.color_picker("Inserts color:", value=col2)

    # To store the buffers that will hold the images
    svg_buffers = []

    for genome in genome_choices:
        st.subheader(genome, divider="gray")
        inserts = comparison[genome].get_by_seq_id(
            seq_ids,
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
                show_mappings=show_mappings,
            )

            # Save the svg figure in a buffer
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format="svg")
            img_buffer.seek(0)
            svg_buffers.append(img_buffer)

            # Display the figure in the webapp
            st.pyplot(fig)
            plt.close()

        for clust_idx, insert_ids in clusters2[genome]:
            clust_inserts = comparison[genome].get_by_insert_id(
                insert_ids, insert_type, filter_threshold
            )
            with st.expander(f"Cluster {clust_idx}:"):
                st.write("\n".join([f"- {x:short}" for x in clust_inserts]))

            # Values chosen based on how well plots scaled in webapp
            h_ratios = (2 + len(clust_inserts), 7)
            figsize = (10, 10 * (h_ratios[0] / 7))

            fig, axs = plt.subplots(2, 1, figsize=figsize, height_ratios=h_ratios)
            fig.suptitle(f"Cluster {clust_idx}")

            axs = comparison[genome].plot_inserts(
                insert_ids=insert_ids,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                buffer=buffer,
                axs=axs,
                feature_types=feature_types,
                colorbar=colorbar,
                col1=col1,
                col2=col2,
            )

            # Save the svg figure in a buffer
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format="svg")
            img_buffer.seek(0)
            svg_buffers.append(img_buffer)

            # Display the figure in the web app
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
        for i, img_buffer in enumerate(svg_buffers):
            zipf.writestr(f"figure_{i+1}.svg", img_buffer.getvalue())

    st.download_button(
        label="Download all figures above",
        key="seq_view_figures_download",
        data=zip_buffer,
        file_name="seq_view.zip",
        mime="image/svg",
        use_container_width=True,
    )


def plot_genomes(
    comparison,
    genome_choice,
    seq_ids,
    clusters,
    insert_type,
    filter_threshold,
    **kwargs,
):
    contig_labels, show_titles = True, True
    num_cols = None

    # Get plotting options
    with st.expander("Plotting options"):
        contig_labels = st.toggle("Show contig/seq labels", value=True)
        cluster_labels = st.toggle("Show cluster labels", value=True)
        show_titles = st.toggle("Show genome labels", value=True)
        facet = st.toggle("Separate plot for each genome", value=False)
        if facet:
            num_cols = st.number_input("Number of columns", value=3, min_value=1)
        palette = st.selectbox("Palette", options=list(color_sequences.keys()))

    if not len(seq_ids):
        seq_ids = None

    # Set seq_labels based on whether some seqs were selected
    seq_labels = defaultdict(dict)
    if contig_labels and (seq_ids is not None or clusters is not None):
        seq_labels = comparison.get_labels(seq_ids, clusters.insert_ids)

    # If cluster labels are preferred, overwrite seq_labels for the sequences
    # in a cluster with the label of a cluster
    if cluster_labels and len(clusters):
        for g, all_labels in clusters.insert_labels.items():
            seq_labels[g].update(all_labels)

    fig = comparison.plot(
        seq_ids={g: seq_ids for g in genome_choice},
        insert_ids=clusters.insert_ids,
        seq_labels=seq_labels,
        show_contig_labels=contig_labels,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        show_titles=show_titles,
        facet_wrap=num_cols,
        palette=palette,
    )

    # Save the svg figure in a buffer
    img = BytesIO()
    fig.savefig(img, format="svg")
    img.seek(0)

    # Show the figure in the web app
    st.pyplot(fig, use_container_width=True)
    plt.close()

    # Show the download button for the svg figure
    st.download_button(
        key="genome_view_figures_download",
        label="Download figure",
        data=img,
        file_name="genome_view.svg",
        mime="image/svg",
        use_container_width=True,
    )

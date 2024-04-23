import argparse

import matplotlib.pyplot as plt

from .main import Pipeline
from .plotting import plot_insert, plot_insert_dists, plot_on_genome


def get_args(*args):
    # Create the CLI using argparse
    # TODO: maybe switch from argparse to click
    parser = argparse.ArgumentParser(description="Generating plots")

    subparsers = parser.add_subparsers(help="sub-command help")

    parser_ncbi = subparsers.add_parser("ncbi", help="Download genome from NCBI")
    parser_ncbi.add_argument("search_term", help="search term for BLAST")
    parser_ncbi.add_argument("email", help="Your email address - needed by BLAST")
    parser_ncbi.add_argument(
        "--retmax",
        default=200,
        help="Maximum number of records to download from NCBI",
    )

    parser_genome = subparsers.add_parser("file", help="Specify a local genome.")
    parser_genome.add_argument("path_to_genome")

    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("output_prefix", help="Output directory")

    parser.add_argument(
        "--fwd_suffix", default="_F", help="Suffix for forward seq. default='_F'"
    )
    parser.add_argument(
        "--rev_suffix", default="_R", help="Suffix for reverse seq. default='_R'"
    )
    parser.add_argument(
        "--insert_types", help="matched, unmatched or both inserts", default="both"
    )
    parser.add_argument(
        "--filter",
        default=0.5,
        help="Threshold value for filtering out blast results with low coverage",
    )

    # Get all the command line arguments
    return parser.parse_args(*args)


def main():
    args = get_args()
    if "path_to_genome" in args:
        search_term = None
        email = None
        retmax = None
        path_to_genome = args.path_to_genome
    else:
        search_term = args.search_term
        email = args.email
        retmax = int(args.retmax)
        path_to_genome = None

    pipeline = Pipeline(
        args.seq_file,
        args.output_prefix,
        genome_file=path_to_genome,
        search_term=search_term,
        email=email,
        retmax=retmax,
        fwd_suffix=args.fwd_suffix,
        rev_suffix=args.rev_suffix,
    )

    # Create linear plots of inserts with annotations
    for seq_id in pipeline.seq_ids:
        inserts = pipeline.get_inserts(
            seq_id, insert_types=args.insert_types, filter_threshold=args.filter
        )
        for i, insert in enumerate(inserts):
            fig, axs = plt.subplots(2, 1, figsize=(10, 10), height_ratios=[3, 5])
            axs = plot_insert(insert, pipeline.genome, axs=axs)
            fig.savefig(pipeline.work_dir / f"{seq_id}_hit{i}.png", dpi=300)
            plt.close()

    # Get all inserts (used in both next plots)
    inserts = pipeline.get_all_inserts(args.insert_types, filter_threshold=args.filter)

    # Create a plot of genome / all contigs as circular plot with inserts
    # layered on top
    fig, ax = plt.subplots(figsize=(10, 30))
    ax = plot_on_genome(inserts, pipeline.genome, ax=ax)
    fig.savefig(pipeline.work_dir / "genome_plot.png", dpi=300)

    # Plots of distributions of insert sizes
    fig, axs = plt.subplots(1, 3, figsize=(12, 5))
    axs = plot_insert_dists(inserts, axs=axs)
    fig.savefig(pipeline.work_dir / "insert_length_dist.png", dpi=300)

    df = pipeline.to_dataframe(
        insert_types=args.insert_types, filter_threshold=args.filter
    )
    df.to_csv(pipeline.work_dir / "inserts.csv", index=False)


if __name__ == "__main__":
    main()

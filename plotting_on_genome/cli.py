import argparse
import matplotlib.pyplot as plt

from .main import Pipeline


def get_args():
    # Create the CLI using argparse
    # TODO: maybe switch from argparse to click
    parser = argparse.ArgumentParser(description="Generating plots")
    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("search_term", help="search term for BLAST")
    parser.add_argument("email", help="Your email address - needed by BLAST")
    parser.add_argument("output_prefix", help="Output directory")
    parser.add_argument(
        "--fwd_suffix", default="_F", help="Suffix for forward seq. default='_F'"
    )
    parser.add_argument(
        "--rev_suffix", default="_R", help="Suffix for reverse seq. default='_R'"
    )
    parser.add_argument(
        "--output", help="matched, unmatched or both inserts", default="both"
    )
    parser.add_argument(
        "--filter",
        default=0.5,
        help="Threshold value for filtering out blast results with low coverage",
    )
    parser.add_argument(
        "--retmax",
        default=200,
        help="Maximum number of records to download from NCBI",
    )

    # Get all the command line arguments
    return parser.parse_args()


def main():
    args = get_args()

    # Create the Pipeline object
    pipeline = Pipeline(
        args.seq_file,
        args.search_term,
        args.email,
        args.output_prefix,
        retmax=int(args.retmax),
        fwd_suffix=args.fwd_suffix,
        rev_suffix=args.rev_suffix,
    )

    # Create linear plots of inserts with annotations
    for seq_id in pipeline.seq_ids:
        for i, insert in enumerate(
            pipeline.get_inserts(
                seq_id, output=args.output, filter_threshold=args.filter
            )
        ):
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            pipeline.plot_insert(insert, axs=axs)
            fig.savefig(pipeline.work_dir / f"{seq_id}_hit{i}.png")
            plt.close()

    # Create a plot of genome / all contigs as circular plot with inserts
    # layered on top
    fig, ax = plt.subplots(figsize=(10, 20))
    pipeline.plot_all_db_seqs(output=args.output, filter_threshold=args.filter, ax=ax)
    fig.savefig(pipeline.work_dir / "genome_plot.png")

    fig, axs = plt.subplots(1, 3, figsize=(12, 5))
    pipeline.plot_insert_dists(
        output=args.output, filter_threshold=args.filter, axs=axs
    )
    fig.savefig(pipeline.work_dir / "insert_length_dist.png")

    pipeline.to_dataframe(output=args.output, filter_threshold=args.filter).to_csv(
        pipeline.work_dir / "inserts.csv", index=False
    )


if __name__ == "__main__":
    main()

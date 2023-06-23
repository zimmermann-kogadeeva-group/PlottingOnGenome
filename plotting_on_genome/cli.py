import argparse
import matplotlib.pyplot as plt

from .main import Pipeline


def main():
    parser = argparse.ArgumentParser(description="Generating plots")
    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("search_term", help="search term for BLAST")
    parser.add_argument("email", help="Your email address - needed by BLAST")
    parser.add_argument("output_prefix", help="Output directory")
    parser.add_argument("--fwd_suffix", default="_F", help="Suffix for forward seq. default='_F'")
    parser.add_argument("--rev_suffix", default="_R", help="Suffix for reverse seq. default='_R'")
    parser.add_argument(
        "--output", help="matched, unmatched or both inserts", default="both"
    )

    args = parser.parse_args()

    pipeline = Pipeline(
        args.seq_file,
        args.search_term,
        args.email,
        args.output_prefix,
        fwd_suffix=args.fwd_suffix,
        rev_suffix=args.rev_suffix
    )

    for seq_id in pipeline.seq_ids:
        for i, insert in enumerate(pipeline.get_inserts(seq_id)):
            fig, axs = plt.subplots(1, 2, figsize=(10, 8))
            pipeline.plot_insert(insert, axs=axs)
            fig.savefig(pipeline.work_dir / f"{seq_id}_hit{i}.png")
            plt.close()

    fig, ax = plt.subplots(figsize=(10, 30))
    pipeline.plot_all_db_seqs(ax=ax)
    fig.savefig(pipeline.work_dir / "genome_plot.png")


if __name__ == "__main__":
    main()

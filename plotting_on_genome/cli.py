import argparse
import matplotlib.pyplot as plt

from .main import Pipeline


def main():
    parser = argparse.ArgumentParser(description="Generating plots")
    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("search_term", help="search term for BLAST")
    parser.add_argument("email", help="Your email address - needed by BLAST")
    parser.add_argument("output_prefix", help="Output directory")
    parser.add_argument("--fwd_suffix", help="Suffix for forward seq.")
    parser.add_argument("--rev_suffix", help="Suffix for reverse seq.")
    parser.add_argument(
        "--output", help="matched, unmatched or both inserts", default="both"
    )

    args = parser.parse_args()

    pipeline = Pipeline(
        args.seq_file,
        args.search_term,
        args.email,
        args.output_prefix,
    )

    for seq_id in pipeline.seq_ids:
        pipeline.plot_all_inserts(seq_id, args.output, save_fmt="png")
        plt.close()

    pipeline.plot_all_db_seqs(args.output, save_fmt="png")

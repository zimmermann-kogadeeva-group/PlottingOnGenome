import argparse
import matplotlib.pyplot as plt

from .main import Pipeline


def main():
    parser = argparse.ArgumentParser(description="Generating plots")
    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("search_term", help="search term for BLAST")
    parser.add_argument("email", help="Your email address - needed by BLAST")
    parser.add_argument("output_prefix", help="Output directory")

    args = parser.parse_args()

    pipeline = Pipeline(
        args.seq_file,
        args.search_term,
        args.email,
        args.output_prefix,
    )
    
    for seq_id in pipeline.blast_results.keys():
        pipeline.plot_all_hsp(seq_id, save_fmt="png")
        plt.close()

    pipeline.plot_all_db_seqs(save_fmt="png")

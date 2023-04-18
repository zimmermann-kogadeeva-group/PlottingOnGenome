
import argparse
from .main import run_pipeline


def main():
    parser = argparse.ArgumentParser(description="Generating plots")
    parser.add_argument("seq_file", help="file with the sequences")
    parser.add_argument("search_term", help="search term for BLAST")
    parser.add_argument("email", help="Your email address - needed by BLAST")
    parser.add_argument("output_prefix", help="Output directory")
    parser.add_argument("--images-prefix", help="(Optional) Folder for images")

    args = parser.parse_args()

    run_pipeline(args.seq_file, args.search_term, args.email, args.output_prefix, args.images_prefix)


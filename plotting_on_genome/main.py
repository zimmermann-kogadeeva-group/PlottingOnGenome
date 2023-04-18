
from pathlib import Path

from .get_genome import main as get_genome
from .run_blast import main as run_blast
from .annotate_sequences import main as annotate_sequences
from .plotting import main as plotting


def run_pipeline(seq_file, search_term, email, output_prefix, images_prefix=None):

    if images_prefix is None:
        images_prefix = output_prefix

    output_prefix = Path(output_prefix)
    if not output_prefix.parent.is_dir():
        print(f"Invalid path: {output_prefix.parent} does not exist!")
        exit()
    if output_prefix.is_dir():
        output_prefix = output_prefix / "output"

    images_prefix = Path(images_prefix)
    if not images_prefix.parent.is_dir():
        print(f"Invalid path: {images_prefix.parent} does not exist!")
        exit()
    if images_prefix.is_dir():
        images_prefix = images_prefix / "output"
    images_prefix = str(images_prefix)

    genome_file = str(output_prefix) + "_combined_contigs.fasta"
    blast_file = str(output_prefix) + "_blast_alignment.txt"
    locus_file = str(output_prefix) + "_combined_locus_tags.csv"
    annot_file = str(output_prefix) + "_annotated.tsv"
    
    get_genome(search_term, email, output_prefix)
    
    run_blast(genome_file, seq_file, blast_file)
    
    annotate_sequences(locus_file, blast_file, annot_file)
    
    plotting(genome_file, annot_file, locus_file, images_prefix)


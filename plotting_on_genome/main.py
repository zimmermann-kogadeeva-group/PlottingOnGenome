
from get_genome import main as get_genome
from run_blast import main as run_blast
from annotate_sequences import main as annotate_sequences
from plotting import main as plotting

def run_pipeline(seq_file, search_term, email, output_prefix, images_prefix=None):

    if images_prefix is None:
        images_prefix = output_prefix

    genome_file = f"{output_prefix}_combined_contigs.fasta"
    blast_file = f"{output_prefix}_blast_alignment.txt"
    locus_file = f"{output_prefix}_combined_locus_tags.csv"
    annot_file = f"{output_prefix}_annotated.tsv"
    
    get_genome(search_term, email, output_prefix)
    
    run_blast(genome_file, seq_file, blast_file)
    
    annotate_sequences(locus_file, blast_file, annot_file)
    
    plotting(genome_file, annot_file, locus_file, images_prefix)


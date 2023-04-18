#!/usr/bin/env python

import subprocess


def main(genome_file, seq_file, output_file):
    # make blast database
    cmd = f"makeblastdb -in {genome_file} -parse_seqids -dbtype nucl"
    subprocess.run(
        cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    # align input sequences with query using blastn with default options
    cmd = f"blastn -query {seq_file} -db {genome_file} -out {output_file} -outfmt 6"
    subprocess.run(
        cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

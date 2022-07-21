#!/usr/bin/env python

import subprocess


def main(genome_file, seq_file, output_file):
    
    # make blast database
    cmd = f"makeblastdb -in {genome_file} -parse_seqids -dbtype nucl"
    subprocess.call(cmd, shell=True) #errorcode = 0 means OK
    
    #align input sequences with query using blastn with default options
    cmd = f"blastn -query {seq_file} -db {genome_file} -out {output_file} -outfmt 6"
    subprocess.call(cmd, shell=True)


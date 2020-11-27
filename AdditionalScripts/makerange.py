#!/usr/bin/env python3
import sys
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""
split a genome into a set of overlapping segments
adapted from nanopolish_makerange.py https://github.com/jts/nanopolish to CulebrONT 
"""
# DRAFT = snakemake.input.draft
# SEGMENT_LENGTH = snakemake.params.segment_len if not snakemake.params.segment_len == '' else 50000
# OVERLAP_LENGTH = snakemake.params.overlap_len if not snakemake.params.overlap_len == '' else 200
# MIN_SEGMENT_LENGTH = 5 * OVERLAP_LENGTH
# OUTPUTFILE = snakemake.output.segments_list

parser = argparse.ArgumentParser(description='Partition a genome into a set of overlapping segments')
parser.add_argument('--draft', help="draft to segment")
parser.add_argument('--segment-length', type=int, default=50000)
parser.add_argument('--overlap-length', type=int, default=200)
parser.add_argument('--output-file', help="output file")
args = parser.parse_args()

DRAFT = args.draft
SEGMENT_LENGTH = args.segment_length
OVERLAP_LENGTH = args.overlap_length
OUTPUT_FILE = args.output_file
MIN_SEGMENT_LENGTH = 5 * OVERLAP_LENGTH

with open(DRAFT, "r") as draft:
    recs = [(title.split(None, 1)[0], len(seq))
            for title, seq in SimpleFastaParser(draft)]

with open(OUTPUT_FILE, "w") as output:
    for name, length in recs:
        n_segments = (length / SEGMENT_LENGTH) + 1
        start = 0
        while start < length:
            end = start + SEGMENT_LENGTH
            # If this segment will end near the end of the contig, extend it to end
            if length - end < MIN_SEGMENT_LENGTH:
                output.write(f"{name}:{start}-{length-1}\n")
                start = length
            else:
                nlen = end + OVERLAP_LENGTH
                output.write(f"{name}:{start}-{nlen}\n")
                start = end
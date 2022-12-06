#!/usr/bin/env python

#-----------------------------------------------------------------------
# data_processing.py
# Author: George Tziampazis, Emmanuel Mhrous
# Takes sequence data and aligns to reference sequence
#-----------------------------------------------------------------------
import sys
import sqlite3
from contextlib import closing
import pandas
from Bio import SeqIO
from Bio import Align
import data_processing


DATABASE_URL = "file:sequence.sqlite?mode=ro"
refseq_file = ["refseq NC_045512.fasta"]
refseq_file = [f"FASTA/{file_name}" for file_name in refseq_file]
refseq_df = data_processing.read_sequences(refseq_file)
refseq_str = refseq_df.sequence[0]

aligner = Align.PairwiseAligner()
aligner.match_score = 1.0
aligner.open_gap_score = -10
aligner.extend_gap_score = -.5

#-----------------------------------------------------------------------

# Returns an alignment object of the reference sequence with the query # sequence
def align_to_refseq(query_seq):
    alignments = aligner.align(refseq_str, query_seq)
    for alignment in alignments:
        return alignment


def main():
    # getting main 
    fasta_file_list = ["part_1_gisaid.fasta", "part_2_gisaid.fasta",\
        "part_3_gisaid.fasta", "part_4_gisaid.fasta",\
            "part_5_gisaid.fasta", "part_6_gisaid.fasta"]
    fasta_file_list = [f"FASTA/{file_name}" for file_name in fasta_file_list]
    sequence_df = data_processing.read_sequences(fasta_file_list)

    aligned = align_to_refseq(sequence_df.sequence[5000])
    print(aligned.substitutions)

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

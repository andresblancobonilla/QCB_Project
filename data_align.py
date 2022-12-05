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
import data_processing


DATABASE_URL = "file:sequence.sqlite?mode=ro"

# new change

#-----------------------------------------------------------------------

# Returns a df with one column of all the sequences 
# in the list of input fasta files
def align_sequences(fasta_file_list):
    # SeqIO.index_db("sequence.idx", fasta_file_list, "fasta")
    sequence_list = []
    
    for fasta_file in fasta_file_list:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            # print(seq_record.seq)
            sequence_list.append(seq_record.seq)
    
    sequence_df = pandas.DataFrame({"sequence":sequence_list})
    return sequence_df


def main():

    refseq_file = ["refseq NC_045512.fasta"]
    refseq_file = [f"FASTA/{file_name}" for file_name in refseq_file]
    sequence_df = data_processing.read_sequences(refseq_file)
    refseq_str = sequence_df.sequence[0]
    print(refseq_str)
    print(sequence_df)
    print(len(sequence_df.index))

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

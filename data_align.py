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
import Bio
from Bio import SeqIO
from Bio import Align
import data_processing
import pickle


DATABASE_URL = "file:sequence.sqlite?mode=ro"
refseq_file = ["refseq NC_045512.fasta"]
refseq_file = [f"FASTA/{file_name}" for file_name in refseq_file]
refseq_df = data_processing.read_sequences(refseq_file)
refseq_str = refseq_df.sequence[0]

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1.0
aligner.open_gap_score = -10
aligner.extend_gap_score = -.5

#-----------------------------------------------------------------------

# Returns an alignment object of the reference sequence with the query # sequence
def align_to_refseq(query_seq):
    alignments = aligner.align(refseq_str, query_seq)
    return alignments[0]

def pickler(data, filepath):
    file = open(filepath, 'wb')
    pickle.dump(data, file)
    file.close()

def depickler(filepath):
    file = open(filepath, 'rb')
    data = pickle.load(file)
    file.close()
    return data

def main():
    # getting main 
    fasta_file_list = ["part_1_gisaid.fasta", "part_2_gisaid.fasta",\
        "part_3_gisaid.fasta", "part_4_gisaid.fasta",\
            "part_5_gisaid.fasta", "part_6_gisaid.fasta"]
    fasta_file_list = [f"FASTA/{file_name}" for file_name in fasta_file_list]
    sequence_df = data_processing.read_sequences(fasta_file_list)
    print(len(sequence_df.index))

    all_alignments = []
    aligned_tuples = []
    for i in range(500):
       aligned = align_to_refseq(sequence_df.sequence[i])
       all_alignments.append(aligned)
       aligned_tuples.append((aligned[0],aligned[1]))

    pickler(all_alignments, "C:/Users/gt512\Documents/Princeton/Code Repos/QCB_Project/pickled_alignments")

    pickler(aligned_tuples, "C:/Users/gt512\Documents/Princeton/Code Repos/QCB_Project/pickled_tuples")

    # depickling file as a variable
    test = depickler("C:/Users/gt512\Documents/Princeton/Code Repos/QCB_Project/pickled_alignments")

    # substitution matrix of first alignment object
    print(test[0].substitutions)



    


#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

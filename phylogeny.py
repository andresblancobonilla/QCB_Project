# -----------------------------------------------------------------------
# mutation_extraction.py
# Author: Fernando Del Castillo, Andres Blanco Bonilla
# Randomly selects 2000 sequences from imported entire sequence file to build phologeny
# Builds a fasta file with 2000 sequences of interest
# Imports aligned output from MAFFT MSA, trims 3' and 5' ends of the virus genome each sample
# Plots some summary statistics of mutations occuring per base pair position, etc.
# with the output from BEAST MCC tree analysis
# -----------------------------------------------------------------------

import sys
from contextlib import closing
import pandas
from Bio import SeqIO
from numpy import random


# Import individual fasta sequences into dataframe with two values: name of sequence (entire string), and the 

def read_sequences(fasta_file_list):
    # SeqIO.index_db("sequence.idx", fasta_file_list, "fasta")
    sequence_list = []
    label = []
    
    for fasta_file in fasta_file_list:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            label.append(seq_record.id)
            sequence_list.append(seq_record.seq)
    
    sequence_df = pandas.DataFrame({"label":label,"sequence":sequence_list})
    return sequence_df

def select_sequences(number,seqdf):
    indices = random.randint(200, size=(number))
    truncatedseqdf = seqdf.iloc[indices]
    return truncatedseqdf

def return_fasta(sequencedf):
    output = ""
    for i in range(1,2000):
        output=output +"\n>"+sequencedf.iloc[i,0]+"\n"+sequencedf.iloc[i,1]
    return output


def main():

    fasta_file_list = ["part_1_gisaid.fasta", "part_2_gisaid.fasta",\
        "part_3_gisaid.fasta", "part_4_gisaid.fasta",\
            "part_5_gisaid.fasta", "part_6_gisaid.fasta"]
    fasta_file_list = [f"FASTA/{file_name}" for file_name in fasta_file_list]
    sequence_df = read_sequences(fasta_file_list)
    # print(sequence_df.iloc[0,:])
    # print(sequence_df)
    truncated_df = select_sequences(2000,sequence_df)
    output_fasta = return_fasta(truncated_df)
    print(output_fasta)
    
    
    
    
    
 

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

#!/usr/bin/env python

#-----------------------------------------------------------------------
# data_processing.py
# Author: Andres Blanco Bonilla
# Reads sequence data into Python for further analysis.
#-----------------------------------------------------------------------
import sys
from Bio import SeqIO

#-----------------------------------------------------------------------

def read_sequences(file_name):
    sequence_list = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        print(seq_record.seq)
        sequence_list.append(seq_record.seq)



def main():
    if len(sys.argv) != 2:
        print('Usage: python data_processing.py fasta_file', file=sys.stderr)
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence_list = read_sequences(fasta_file)
    print(sequence_list)

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

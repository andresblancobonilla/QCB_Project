#!/usr/bin/env python

#-----------------------------------------------------------------------
# data_processing.py
# Author: Andres Blanco Bonilla
# Reads sequence data into Python for further analysis.
#-----------------------------------------------------------------------
import sys
import sqlite3
from contextlib import closing
import pandas
from Bio import SeqIO


DATABASE_URL = "file:sequence.sqlite?mode=ro"

#-----------------------------------------------------------------------

# Returns a df with one column of all the sequences 
# in the list of input fasta files
def read_sequences(fasta_file_list):
    # SeqIO.index_db("sequence.idx", fasta_file_list, "fasta")
    sequence_list = []
    
    for fasta_file in fasta_file_list:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            # print(seq_record.seq)
            sequence_list.append(seq_record.seq)
    
    sequence_df = pandas.DataFrame({"sequence":sequence_list})
    return sequence_df

def get_location_id_date():
    try:
        with sqlite3.connect(DATABASE_URL, isolation_level=None,
                             uri=True) as connection:
            with closing(connection.cursor()) as cursor:
                query = "SELECT key FROM offset_data;"
                cursor.execute(query)
                id_list = []
                row = cursor.fetchone()
                while row is not None:
                    id_list.append(row[0])
                return id_list

    except Exception as ex:
        print(sys.argv[0] + ": " + str(ex), file=sys.stderr)
        sys.exit(1)



def main():

    fasta_file_list = ["part_1_gisaid.fasta", "part_2_gisaid.fasta",\
        "part_3_gisaid.fasta", "part_4_gisaid.fasta",\
            "part_5_gisaid.fasta", "part_6_gisaid.fasta"]
    fasta_file_list = [f"FASTA/{file_name}" for file_name in fasta_file_list]
    sequence_df = read_sequences(fasta_file_list)
    print(sequence_df)
    print(len(sequence_df.index))
    # print(list(sequence_df["sequence"])[3])

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

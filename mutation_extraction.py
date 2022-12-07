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
import data_align
import pickle
from itertools import islice
import collections

#-----------------------------------------------------------------------

GENOME_SEGMENTS = {"ORF1ab": (265,21554), "S":(21562,25383), "ORF3a":(25392,26219),"E":(26244,26471),\
    "M":(26522,27190), "ORF6":(27201,27386), "ORF7a":(27393,27758), "ORF7b":(27755,27886), "ORF8":(27893,28258),\
        "ORF10":(29557,29673), "N":(28273,29532)}

CODON_TABLE = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G',
    '---': 'none' # None
}



#-----------------------------------------------------------------------

def consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)

#-----------------------------------------------------------------------

def count_mutations():
    aligned_sequences_list = data_align.depickler("pickled_tuples")
    base_pair_changes_dict = {}
    base_pair_changes_position_dict = {}
    segment_mutations_dict = {"ORF1ab": {"mis":0, "syn":0}, "S":{"mis":0, "syn":0}, "ORF3a":{"mis":0, "syn":0},\
        "E":{"mis":0, "syn":0},"M":{"mis":0, "syn":0}, "ORF6":{"mis":0, "syn":0},\
        "ORF7":{"mis":0, "syn":0}, "ORF8":{"mis":0, "syn":0}, "ORF10":{"mis":0, "syn":0}, "N":{"mis":0, "syn":0}}

    for reference, query in aligned_sequences_list:
        print(len(reference))
        print(reference[264:270])
        print(reference[21551:21557])
        notprinted = True
        indices = iter(range(len(reference)))
        for index in indices:
            current_segment = "Intergenic"
            for segment, segment_range in GENOME_SEGMENTS.items():
                segment_start = segment_range[0]
                segment_end = segment_range[1]

                if index == segment_range[0]:
                    current_segment = segment
                    current_ref_segment = reference[segment_range[0]:segment_range[1]]
                    current_query_segment = query[segment_range[0]:segment_range[1]]
                    
                    ref_codons = [current_ref_segment[i:i + 3] for i in range(0, len(current_ref_segment), 3)]
                    query_codons = [current_query_segment[i:i + 3] for i in range(0, len(current_query_segment), 3)]
                    if notprinted:
                        print()
                        print(ref_codons)
                        print(len(ref_codons))
                        print()
                        notprinted = False
                    current_index = index
                    for ref_codon, query_codon in zip(ref_codons, query_codons):
                        for ref_base, query_base in zip(ref_codon, query_codon):
                            if ref_base != query_base:
                                base_change = ref_base + ">" + query_base
                                base_change_position = str(index) + base_change
                                if base_pair_changes_dict.get(base_change):
                                    base_pair_changes_dict[base_change]+=1
                                else:
                                    base_pair_changes_dict[base_change] = 1

                                if base_pair_changes_position_dict.get(base_change_position):
                                    base_pair_changes_position_dict[base_change_position]+=1
                                else:
                                    base_pair_changes_position_dict[base_change_position] = 1

                        if ref_codon != query_codon:
                            if CODON_TABLE[ref_codon] != CODON_TABLE[query_codon]:
                                segment_mutations_dict[current_segment]["mis"]+=1
                            else:
                                segment_mutations_dict[current_segment]["syn"]+=1
                        break
    
                else:
                    break
                        
                        
            if current_segment != "Intergenic":
                consume(indices, segment_end - segment_start)

            else:
                ref_base = reference[index]
                query_base = query[index]
                if ref_base != query_base:
                    base_change = ref_base + ">" + query_base
                    base_change_position = str(index) + base_change

                    if base_pair_changes_dict.get(base_change):
                        base_pair_changes_dict[base_change]+=1
                    else:
                        base_pair_changes_dict[base_change] = 1

                    if base_pair_changes_position_dict.get(base_change_position):
                        base_pair_changes_position_dict[base_change_position]+=1
                    else:
                        base_pair_changes_position_dict[base_change_position] = 1

    return base_pair_changes_dict, base_pair_changes_position_dict, segment_mutations_dict
                
def main():
    changes_tuple = count_mutations()
    data_align.pickler(changes_tuple, "basepairchanges")
    print(changes_tuple[0])
    print()
    print(changes_tuple[1])
    print()
    print(changes_tuple[2])

#-----------------------------------------------------------------------
if __name__ == '__main__':
    main()

#!/usr/bin/env python

# -----------------------------------------------------------------------
# data_processing.py
# Author: George Tziampazis, Emmanuel Mhrous
# Takes sequence data and aligns to reference sequence
# -----------------------------------------------------------------------
import sys
import sqlite3
from contextlib import closing
import pandas
import Bio
from Bio import SeqIO
from Bio import Align
from Bio.Seq import translate
import data_processing
import data_align
import pickle
from itertools import islice
import collections

# -----------------------------------------------------------------------



NON_CODING_SEGMENTS = {"3'-UTR": (29674, 29903), "5'-UTR": (0, 265)}

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
    '---': 'none'  # None
}


# -----------------------------------------------------------------------


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

def adjust_positions(reference, orfs):
    new_orfs = orfs
    start_inserts = 0
    total_inserts = 0
    if "-" in reference[0:265]:
        start_inserts = reference[0:265].count("-")

    total_inserts+=start_inserts
    # reference = reference[265:]
    orf1ab_start_index = orfs["ORF1ab"][0] + total_inserts

    orf1ab_join_index = reference.index("TAAACGG") + 4
    
    orf1ab_end_index = orf1ab_end_index = orfs["ORF1ab"][2] + total_inserts
    
    orf1ab = reference[orf1ab_start_index:orf1ab_end_index]
    # if reference[orf1ab_join_index] != "C":
    #     print("join:" + str(reference[orf1ab_join_index]) + "\n")
    orf1ab_inserts = 0
    if "-" in orf1ab:
        orf1ab_inserts = orf1ab.count("-")
    total_inserts+=orf1ab_inserts
    orf1ab_end_index+=orf1ab_inserts
    
    new_orfs["ORF1ab"] = (orf1ab_start_index, orf1ab_join_index, orf1ab_end_index)

    # "S": (21562, 25384), "ORF3a": (25392, 26220)
    # s_start_index = orfs["S"][0] + total_inserts
    # s_end_index = orfs["S"][1] + total_inserts
    
    # s = reference[s_start_index:s_end_index]
    # s_inserts = 0
    # if "-" in s:
    #     s_inserts = s.count("-")
    # s_end_index+=s_inserts
    # total_inserts+=s_inserts

    # orfs["S"] = (s_start_index, s_end_index)
    
    for orf in list(orfs.keys()):
        if orf == "ORF1ab":
            continue
        orf_start_index = orfs[orf][0] + total_inserts
        orf_end_index = orfs[orf][1] + total_inserts
        
        orf_sequence = reference[orf_start_index:orf_end_index]
        orf_inserts = 0
        if "-" in orf_sequence:
            orf_inserts = orf_sequence.count("-")
        orf_end_index+=orf_inserts
        total_inserts+=orf_inserts
        new_orfs[orf] = (orf_start_index, orf_end_index)
    
    return new_orfs

    # orf1ab_start_index = reference.index("ATG")
    # total_inserts = 0
    # UTR_3_inserts = 0
    # if orf1ab_start_index != 265:
    #     UTR_3_inserts = orf1ab_start_index - 265
    # orf1ab = reference[orf1ab_start_index:]
    # orf1ab_join_index = orf1ab.index("TAAACGG") + 4 + UTR_3_inserts
    # orf1ab_inserts = 0
    # # if orf1ab_join_index != orf1ab_start_index + 13202:
    # #     orf1ab_inserts = orf1ab_join_index - 13467
    # orf1ab_end_index = orf1ab.index("TAA") + UTR_3_inserts
    # orfs["ORF1AB"] = (orf1ab_start_index, orf1ab_end_index, orf1ab_join_index)
#-----------------------------------------------------------------------


def split_orf1ab(full_sequence, range):
    orf1ab_sequence = full_sequence[range[0]:range[1] + 1] + \
        full_sequence[range[1]:range[2]]
    
    # orf1ab_sequence = orf1ab_sequence.replace("-", "N")

    # # print(orf1ab_sequence)
    # try:
    #     translated = translate(orf1ab_sequence)
    # except Exception as ex:
    #     print(ex)
    #     # print(len(orf1ab_sequence))
    #     # print(orf1ab_sequence)
    #     sys.exit(1)

    # print(translated)
    return orf1ab_sequence

# -----------------------------------------------------------------------


def count_mutations():
    aligned_sequences_list = data_align.depickler("pickled_tuples")
    # insertion_count = 0
    # deletion_count = 0
    # notprinted = True
    # for rs, qs in aligned_sequences_list:
    #     if "-" in rs:
    #         if notprinted:
    #             print(rs[-60:])
    #             print(qs[-60:])
    #             print()
    #         insertion_count+=1
    #     if "-" in qs:
    #         deletion_count+=1
    base_pair_changes_dict = {}
    base_pair_changes_position_dict = {}
    # notprinted = True
    sequence_count = 0
    segment_mutations_dict = {"3'-UTR": {"mut": 0, "in": 0, "del": 0}, "5'-UTR": {"mut": 0, "in": 0, "del": 0}, "Intergenic": {"mut": 0, "in": 0, "del": 0}, "ORF1ab": {"mis": 0, "syn": 0, "in": 0, "del": 0},
                              "S": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "ORF3a": {"mis": 0, "syn": 0, "in": 0, "del": 0},
                              "E": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "M": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "ORF6": {"mis": 0, "syn": 0, "in": 0, "del": 0},
                              "ORF7a": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "ORF7b": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "ORF8": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "ORF10": {"mis": 0, "syn": 0, "in": 0, "del": 0}, "N": {"mis": 0, "syn": 0, "in": 0, "del": 0}}
    sequence_index = 0
    sequences_base_pair_changes_position_list = []
    for reference, query in aligned_sequences_list:
        sequence_changes = {"sequence_index":sequence_index, "mutations":[]}
        sequence_index+=1
        end_UTR_5 = 265
        start_UTR_3 = 29674
        orfs = {"ORF1ab": (265, 21555, 13467), "S": (21562, 25384), "ORF3a": (25392, 26220), "E": (26244, 26472),
                   "M": (26522, 27191), "ORF6": (27201, 27387), "ORF7a": (27393, 27759), "ORF7b": (27755, 27887), "ORF8": (27893, 28259),
                   "ORF10": (29557, 29674), "N": (28273, 29533)}
        if "-" in reference:
            orfs = adjust_positions(reference, orfs)
            end_UTR_5 = orfs["ORF1ab"][0]
            start_UTR_3 = orfs["ORF10"][1]
            if reference[start_UTR_3 - 3:start_UTR_3+1] != "TAGC":
                print(reference)
                print(orfs)
            # print(reference[start_UTR_3 - 3:start_UTR_3+1])
            

        sequence_count += 1

        # if "":
            # handle_deletions(reference, query)

        indices = iter(range(len(reference)))
        for index in indices:
            current_index = index
            current_segment = "Intergenic"
            for segment, segment_range in orfs.items():
                segment_start = segment_range[0]
                segment_end = segment_range[1]

                if index == segment_start:
                    current_segment = segment
                    if segment == "ORF1ab":
                        current_ref_segment = split_orf1ab(reference, segment_range)
                        current_query_segment = split_orf1ab(query, segment_range)
                    else:
                        current_ref_segment = reference[segment_start:segment_end]
                        current_query_segment = query[segment_start:segment_end]

                    ref_codons = [current_ref_segment[i:i + 3]
                                  for i in range(0, len(current_ref_segment), 3)]
                    query_codons = [current_query_segment[i:i + 3]
                                    for i in range(0, len(current_query_segment), 3)]

                    for ref_codon, query_codon in zip(ref_codons, query_codons):
                        codon_indices = iter(range(len(ref_codon)))
                        if "-" in query_codon or "-" in ref_codon:
                            break
                        for codon_index in codon_indices:
                            ref_base = ref_codon[codon_index]
                            query_base = query_codon[codon_index]
                            if ref_base != query_base:
                                # if query_base == "-":
                                #     break
                                # # if index == 1604:
                                # #     print()
                                # #     print("stop here!")
                                # #     print()
                                # start_index = codon_index
                                # start_current_index = current_index
                                # end_index = codon_index
                                # current_char = "-"
                                # while current_char == "-" and end_index != 2:
                                #     end_index+=1
                                #     current_index+=1
                                #     current_char = query_codon[end_index]

                                # base_change_position = f"{start_current_index + 1}_{current_index + 1}del"
                                # if current_index == 1605:
                                #     print()
                                #     print("here")
                                #     print()

                                # if base_pair_changes_position_dict.get(base_change_position):
                                #     base_pair_changes_position_dict[base_change_position]+=1
                                # else:
                                #     base_pair_changes_position_dict[base_change_position] = 1
                                #     segment_mutations_dict[current_segment]["del"]+=1
                                # if start_index != end_index:
                                #     consume(codon_indices, end_index - start_index)
                                # else:
                                base_change = ref_base + ">" + query_base
                                base_change_position = str(
                                    current_index + 1) + base_change
                                sequence_changes["mutations"].append(base_change_position)

                                if base_pair_changes_position_dict.get(base_change_position):
                                    base_pair_changes_position_dict[base_change_position] += 1
                                else:
                                    base_pair_changes_position_dict[base_change_position] = 1
                                    if base_pair_changes_dict.get(base_change):
                                        base_pair_changes_dict[base_change] += 1
                                    else:
                                        base_pair_changes_dict[base_change] = 1

                                    if ref_codon != query_codon:
                                        if not CODON_TABLE.get(query_codon) or CODON_TABLE[ref_codon] != CODON_TABLE[query_codon]:
                                            segment_mutations_dict[current_segment]["mis"] += 1
                                        else:
                                            segment_mutations_dict[current_segment]["syn"] += 1
                            current_index += 1
                    break

                else:
                    continue

            if current_index != index:
                consume(indices, current_index - index)

            if current_index < end_UTR_5:
                current_segment = "5'-UTR"
            elif current_index >= start_UTR_3:
                current_segment = "3'-UTR"

            ref_base = reference[current_index]
            query_base = query[current_index]
            if ref_base != query_base:
                if query_base == "-":
                    start_index = current_index
                    end_index = current_index
                    current_char = "-"
                    while current_char == "-" and end_index != len(query) - 1:
                        end_index += 1
                        current_char = query[end_index]
                    if start_index != end_index:
                        base_change_position = f"{start_index+1}_{end_index}del"
                    else:
                        base_change_position = f"{start_index+1}del"

                    if base_pair_changes_position_dict.get(base_change_position):
                        base_pair_changes_position_dict[base_change_position] += 1
                    else:
                        base_pair_changes_position_dict[base_change_position] = 1
                        segment_mutations_dict[current_segment]["del"] += 1
                    if start_index != end_index:
                        consume(indices, end_index -
                                start_index - 1)
                elif ref_base == "-":
                    start_index = current_index
                    end_index = current_index
                    current_char = "-"
                    while current_char == "-" and end_index != len(reference) - 1:
                        end_index += 1
                        current_char = reference[end_index]
                    if start_index != end_index:
                        base_change_position = f"{start_index+1}_{end_index}in"
                    else:
                        base_change_position = f"{start_index+1}in"
                        

                    if base_pair_changes_position_dict.get(base_change_position):
                        base_pair_changes_position_dict[base_change_position] += 1
                    else:
                        base_pair_changes_position_dict[base_change_position] = 1
                        segment_mutations_dict[current_segment]["in"] += 1
                    if start_index != end_index:
                        consume(indices, end_index -
                                start_index - 1)
                else:
                    base_change = ref_base + ">" + query_base
                    base_change_position = str(
                        index + 1) + base_change

                    if base_pair_changes_dict.get(base_change):
                        base_pair_changes_dict[base_change] += 1
                    else:
                        base_pair_changes_dict[base_change] = 1

                    if base_pair_changes_position_dict.get(base_change_position):
                        base_pair_changes_position_dict[base_change_position] += 1
                    else:
                        base_pair_changes_position_dict[base_change_position] = 1
                        segment_mutations_dict[current_segment]["mut"] += 1
                sequence_changes["mutations"].append(base_change_position)
        sequences_base_pair_changes_position_list.append(sequence_changes)

    print(sequence_count)
    print()
    
    
    sorted_base_pair_changes = sorted(base_pair_changes_dict.items(), key=lambda x:x[1], reverse=True)
    base_pair_changes_dict = dict(sorted_base_pair_changes)

    sorted_base_pair_changes_position = sorted(base_pair_changes_position_dict.items(), key=lambda x:x[1], reverse=True)
    base_pair_changes_position_dict = dict(sorted_base_pair_changes_position)


    return base_pair_changes_dict, base_pair_changes_position_dict, segment_mutations_dict, sequences_base_pair_changes_position_list


def main():
    changes_tuple = count_mutations()
    data_align.pickler(changes_tuple, "basepairchanges")
    print(changes_tuple[0])
    print()
    print(changes_tuple[1])
    print()
    print(changes_tuple[2])
    print()
    print(changes_tuple[3])


# -----------------------------------------------------------------------
if __name__ == '__main__':
    main()

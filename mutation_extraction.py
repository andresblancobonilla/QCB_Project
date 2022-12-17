#!/usr/bin/env python

#-----------------------------------------------------------------------
# mutation_extraction.py
# Author: Andres Blanco Bonilla, Fernando Del Castillo
# Extracts mutations out of sequence by comparing it to reference
#-----------------------------------------------------------------------

import Bio
from Bio import SeqIO
from Bio import Align
from Bio.Seq import translate
import data_processing
import data_align
import pickle
from itertools import islice
import collections

#-----------------------------------------------------------------------
# This function is copied from the itertools Python documentation.
# https://docs.python.org/3/library/itertools.html#recipes
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

# If an insertion is encountered sequence,
# this function adjusts the windows of the ORFS as needed.
def adjust_positions(reference, orfs):
    new_orfs = orfs
    start_inserts = 0
    total_inserts = 0
    if "-" in reference[0:265]:
        start_inserts = reference[0:265].count("-")

    total_inserts += start_inserts
    orf1ab_start_index = orfs["ORF1ab"][0] + total_inserts
    
    # This accounts for the -1 frameshift
    # at sequence position 13468 (here index 13467).
    # The "C" nucleotide is essentially read twice.
    orf1ab_join_index = reference.index("TAAACGG") + 4

    orf1ab_end_index = orf1ab_end_index = orfs["ORF1ab"][2] + \
        total_inserts

    orf1ab = reference[orf1ab_start_index:orf1ab_end_index]
    orf1ab_inserts = 0
    if "-" in orf1ab:
        orf1ab_inserts = orf1ab.count("-")
    total_inserts += orf1ab_inserts
    orf1ab_end_index += orf1ab_inserts

    new_orfs["ORF1ab"] = (orf1ab_start_index,
                          orf1ab_join_index, orf1ab_end_index)

    for orf in list(orfs.keys()):
        if orf == "ORF1ab":
            continue
        orf_start_index = orfs[orf][0] + total_inserts
        orf_end_index = orfs[orf][1] + total_inserts

        orf_sequence = reference[orf_start_index:orf_end_index]
        orf_inserts = 0
        if "-" in orf_sequence:
            orf_inserts = orf_sequence.count("-")
        orf_end_index += orf_inserts
        total_inserts += orf_inserts
        new_orfs[orf] = (orf_start_index, orf_end_index)

    return new_orfs

#-----------------------------------------------------------------------
# This accounts for the -1 frameshift that occurs during the translation
# of ORF1ab. ORF1ab's position is join(266..13468,13468..21555), as per
# the NCBI page for the reference sequence.
def split_orf1ab(full_sequence, range):
    orf1ab_sequence = full_sequence[range[0]:range[2] + 1] + \
        full_sequence[range[2]:range[1]]
    return orf1ab_sequence

#-----------------------------------------------------------------------
# Count both nucleotide and amino acid mutations in every query sequence
# compared to the reference. Generally, check if a mutation, insertion,
# or deletion has already been encountered in a previous sequence
# before considering it as a new variant and incrementing counts.

def count_mutations():
    # Read the pairs of aligned sequences in from our file
    aligned_sequences_list = data_align.depickler("all_pickled_tuples")
    
    base_pair_changes_dict = {}
    base_pair_changes_position_dict = {}
    amino_acid_changes_position_dict = {}
    segment_mutations_dict = {"3'-UTR": {"mut": 0, "in": 0, "del": 0},
                              "5'-UTR": {"mut": 0, "in": 0, "del": 0},
                              "Intergenic": {"mut": 0, "in": 0, "del": 0},
                              "ORF1ab": {"mis": 0, "syn": 0, "in": 0, "del": 0, "ifi": 0, "ifd": 0, "fsd":0, "sg": 0},
                              "S": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "ORF3a": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "E": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "M": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "ORF6": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "ORF7a": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "ORF8": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "ORF10": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0},
                              "N": {"mis": 0, "syn": 0, "in": 0, "del": 0,  "ifi": 0, "ifd": 0, "fsd":0,  "sg": 0}}
    sequences_base_pair_changes_position_list = []
    sequences_amino_acid_changes_position_list = []

    sequence_index = 0
    sequence_count = 0
    
    for reference, query in aligned_sequences_list:
        sequence_changes = {
            "sequence_index": sequence_index, "mutations": []}
        sequence_amino_acid_changes = {
            "sequence_index": sequence_index, "mutations": []}
        sequence_index += 1

        end_UTR_5 = 265
        start_UTR_3 = 29674
        orfs = {"ORF1ab": (265, 21555, 13467), "S": (21562, 25384), "ORF3a": (25392, 26220),
                "E": (26244, 26472), "M": (26522, 27191), "ORF6": (27201, 27387),
                "ORF7a": (27393, 27759), "ORF8": (27893, 28259),
                "N": (28273, 29533), "ORF10": (29557, 29674)}
        
        # Adjust these ORF ranges if there is an insertion
        if "-" in reference:
            orfs = adjust_positions(reference, orfs)
            end_UTR_5 = orfs["ORF1ab"][0]
            start_UTR_3 = orfs["ORF10"][1]

        sequence_count += 1


        indices = iter(range(len(reference)))
        for index in indices:
            current_index = index
            current_segment = "Intergenic"

            for segment, segment_range in orfs.items():
                segment_start = segment_range[0]
                segment_end = segment_range[1]

                if index >= segment_start and index < segment_end:
                    current_segment = segment
                    
                    # if at the start of an ORF, loop through the codons
                    # in the ORF
                    if index == segment_start:
                        if current_segment == "ORF1ab":
                            current_ref_segment = split_orf1ab(
                                reference, segment_range)
                            current_query_segment = split_orf1ab(
                                query, segment_range)
                            # index must be corrected later
                            # to account for the -1 frameshift
                            corrected_index = False
                        else:
                            current_ref_segment = reference[segment_start:segment_end]
                            current_query_segment = query[segment_start:segment_end]

                        ref_codons = [current_ref_segment[i:i + 3]
                                      for i in range(0, len(current_ref_segment), 3)]
                        query_codons = [current_query_segment[i:i + 3]
                                        for i in range(0, len(current_query_segment), 3)]

                        # BioPython translate will fail if gaps ("-")
                        # are in the sequence

                        ref_amino_acids = translate(
                            current_ref_segment.replace("-", ""))
                        query_amino_acids = translate(
                            current_query_segment.replace("-", ""))
                        segment_amino_acid_index = 0

                        for ref_codon, query_codon, ref_amino_acid,\
                        query_amino_acid in zip(ref_codons,\
                        query_codons, ref_amino_acids,\
                        query_amino_acids):

                            segment_amino_acid_index += 1
                            
                            # hanle deletions
                            if "-" in query_codon:

                                if ref_amino_acid != query_amino_acid:   
                                    amino_acid_change_position = f"{ref_amino_acid}{segment_amino_acid_index}del"
                                else:
                                    amino_acid_change_position = f"{ref_amino_acids[segment_amino_acid_index]}{segment_amino_acid_index + 1}del"
                                    
                                sequence_amino_acid_changes["mutations"].append(amino_acid_change_position)
                                    
                                if amino_acid_changes_position_dict.get(amino_acid_change_position):
                                    amino_acid_changes_position_dict[amino_acid_change_position] += 1
                                else:
                                    amino_acid_changes_position_dict[amino_acid_change_position] = 1
                                    if len(current_query_segment) % 3 != 0:
                                        segment_mutations_dict[current_segment]["fsd"] += 1
                                    else:
                                        segment_mutations_dict[current_segment]["ifd"] += 1
                
                                # quit codon loop if insertion is encountered
                                break
                                        
                            # handle insertions
                            if "-" in ref_codon:
                                if ref_amino_acid != query_amino_acid:   
                                    amino_acid_change_position = f"{query_amino_acid}{segment_amino_acid_index}in"
                                else:
                                    amino_acid_change_position = f"{query_amino_acids[segment_amino_acid_index]}{segment_amino_acid_index + 1}in"
                                    
                                sequence_amino_acid_changes["mutations"].append(amino_acid_change_position)
                                                            
                                if amino_acid_changes_position_dict.get(amino_acid_change_position):
                                    amino_acid_changes_position_dict[amino_acid_change_position] += 1
                                else:
                                    amino_acid_changes_position_dict[amino_acid_change_position] = 1
                                    # segment_mutations_dict[current_segment]["mis"] += 1
                                    if len(current_ref_segment) % 3 != 0:
                                        segment_mutations_dict[current_segment]["fsd"] += 1
                                    else:
                                        segment_mutations_dict[current_segment]["ifi"] += 1
                                # quit codon loop if insertion is encountered
                                break
                            
                            # handle stop-gained mutations
                            if ref_amino_acid != query_amino_acid and query_amino_acid == "*":
                                amino_acid_change_position = f"{ref_amino_acid}{segment_amino_acid_index}{query_amino_acid}"
                                    
                                sequence_amino_acid_changes["mutations"].append(amino_acid_change_position)
                                
                                if amino_acid_changes_position_dict.get(amino_acid_change_position):
                                    amino_acid_changes_position_dict[amino_acid_change_position] += 1
                                else:
                                    amino_acid_changes_position_dict[amino_acid_change_position] = 1
                                    segment_mutations_dict[current_segment]["mis"] += 1
                                    segment_mutations_dict[current_segment]["sg"] += 1
                                break
                            
                            # iterate through bases in codons
                            for ref_base, query_base in zip(ref_codon, query_codon):
                                current_index += 1
                                # ignore ambigious reads
                                if ref_base != query_base and query_base in ["A", "C", "T", "G"]:
                                    # correct index to account for ORF1ab frameshit
                                    if current_segment == "ORF1ab"\
                                    and current_index > orfs["ORF1ab"][2]\
                                    and not corrected_index:
                                        current_index-=1
                                        corrected_index = True

                                    base_change = ref_base + ">" + query_base

                                    base_change_position = str(
                                        current_index) + base_change
                                    sequence_changes["mutations"].append(
                                        base_change_position)
                                    
                                    if base_pair_changes_position_dict.get(base_change_position):
                                        base_pair_changes_position_dict[base_change_position] += 1
                                    else:
                                        base_pair_changes_position_dict[base_change_position] = 1
                                        if base_pair_changes_dict.get(base_change):
                                            base_pair_changes_dict[base_change] += 1
                                        else:
                                            base_pair_changes_dict[base_change] = 1

                                    # ignore ambiguous reads
                                    if query_amino_acid != 'X':
                                        if ref_amino_acid != query_amino_acid:
                                                amino_acid_change_position = f"{ref_amino_acid}{segment_amino_acid_index}{query_amino_acid}"
                                                sequence_amino_acid_changes["mutations"].append(
                                                amino_acid_change_position)
                                                if amino_acid_changes_position_dict.get(amino_acid_change_position):
                                                    amino_acid_changes_position_dict[amino_acid_change_position] += 1
                                                else:
                                                    amino_acid_changes_position_dict[amino_acid_change_position] = 1
                                                    segment_mutations_dict[current_segment]["mis"] += 1
                        
                                        else:
                                            if base_pair_changes_position_dict.get(base_change_position) == 1:
                                                segment_mutations_dict[current_segment]["syn"] += 1

                    # exit the segment loop if inside a segment
                    break

                else:
                    continue

            if current_index < end_UTR_5:
                current_segment = "5'-UTR"
            elif current_index >= start_UTR_3:
                current_segment = "3'-UTR"
            
            end_of_sequence_index = len(reference) - 1

            ref_base = reference[current_index]
            query_base = query[current_index]

            # count mutations/variations in non-coding regions
            if ref_base != query_base and query_base in ["A", "C", "T", "G", "-"]:
                
                # handle deletions by iterating through all the deletions
                # that are in a row and then skipping forward in the
                # outer index loop
                if query_base == "-":
                    start_index = current_index
                    end_index = current_index
                    next_char = ""
                    if end_index < end_of_sequence_index:
                        next_char = query[end_index + 1]
                    
                    while next_char == "-" and end_index < end_of_sequence_index:
                        end_index +=1
                        next_char = query[end_index]

                    if start_index != end_index:
                        base_change_position = f"{start_index+1}_{end_index}del"
                    else:
                        base_change_position = f"{start_index+1}del"
                        
                    
                    # ignore deletions if they start at the start of the
                    # sequence or end at the end of the sequence
                    if end_index != end_of_sequence_index and start_index != 0:
                        if base_pair_changes_position_dict.get(base_change_position):
                            base_pair_changes_position_dict[base_change_position] += 1
                        else:
                            base_pair_changes_position_dict[base_change_position] = 1
                            segment_mutations_dict[current_segment]["del"] += 1
                    if start_index != end_index:
                        consume(indices, end_index - start_index - 1)

                # handle insertions by iterating through all the insertions
                # that are in a row and then skipping forward in the
                # outer index loop
                elif ref_base == "-":
                    start_index = current_index
                    end_index = current_index
                    next_char = ""
                    if end_index < end_of_sequence_index:
                        next_char = reference[end_index + 1]
                    
                    while next_char == "-" and end_index < end_of_sequence_index:
                        end_index +=1
                        next_char = reference[end_index]

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
                        consume(indices, end_index - start_index - 1)
                
                # count mutations in non-coding regions only 
                elif current_segment in ["5'-UTR", "3'-UTR", "Intergenic"]:
                    base_change = ref_base + ">" + query_base
                    base_change_position = str(
                        index + 1) + base_change


                    if base_pair_changes_position_dict.get(base_change_position):
                        base_pair_changes_position_dict[base_change_position] += 1
                    else:
                        base_pair_changes_position_dict[base_change_position] = 1
                        segment_mutations_dict[current_segment]["mut"] += 1
                        if base_pair_changes_dict.get(base_change):
                            base_pair_changes_dict[base_change] += 1
                        else:
                            base_pair_changes_dict[base_change] = 1

                sequence_changes["mutations"].append(
                    base_change_position)
            else:
                continue

        sequences_base_pair_changes_position_list.append(
            sequence_changes)

        sequences_amino_acid_changes_position_list.append(
            sequence_amino_acid_changes)

    sorted_base_pair_changes = sorted(
        base_pair_changes_dict.items(), key=lambda x: x[1], reverse=True)
    base_pair_changes_dict = dict(sorted_base_pair_changes)

    sorted_base_pair_changes_position = sorted(
        base_pair_changes_position_dict.items(), key=lambda x: x[1], reverse=True)
    base_pair_changes_position_dict = dict(
        sorted_base_pair_changes_position)

    sorted_amino_acid_changes_position = sorted(
        amino_acid_changes_position_dict.items(), key=lambda x: x[1], reverse=True)
    amino_acid_changes_position_dict = dict(
        sorted_amino_acid_changes_position)

    return base_pair_changes_dict, base_pair_changes_position_dict,\
        segment_mutations_dict, amino_acid_changes_position_dict,\
            sequences_base_pair_changes_position_list,\
            sequences_amino_acid_changes_position_list


def main():
    changes_tuple = count_mutations()
    data_align.pickler(changes_tuple, "all_base_pair_changes_updated_again")
    print(changes_tuple[0])
    print()
    print(changes_tuple[1])
    print()
    print(changes_tuple[2])
    print()
    print(changes_tuple[3])
    print()
    # print(changes_tuple[4])
    # print()
    # print(changes_tuple[5])
    # print()

# -----------------------------------------------------------------------
if __name__ == '__main__':
    main()

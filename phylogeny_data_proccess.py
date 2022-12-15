# -----------------------------------------------------------------------
# mutation_extraction.py
# Author: Fernando Del Castillo, Andres Blanco Bonilla
# Randomly selects 2000 sequences from imported entire sequence file to build phologeny
# Builds a fasta file with 2000 sequences of interest
# Imports aligned output from MAFFT MSA, trims 3' and 5' ends of the virus genome each sample
# Plots some summary statistics of mutations occuring per base pair position, etc.
# with the output from BEAST MCC tree analysis
# -----------------------------------------------------------------------

from Bio import SeqIO

input_file = "mafft_alignment_aln.fasta"

#!/usr/bin/python3

from os import sys
import PyDawg as pd
from PyDawg import Segment
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

trickString = \
"""# Example: Simulate DNA evolution along a tree
# See readme.txt for an explanation on the structure of an input file.
# Simulation results are sent to stdout.

## Tree Section ################################################################
# 
# Use a constant tree.

[Tree]
Tree = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);"

## Subst Section ###############################################################
#
# Use an HKY substitution model with a transition rate of 2.0 and a transversion
# rate of 1.0.  Allele frequences are 0.3 A, 0.2 C, 0.2 G, and 0.3 T.

[Subst]
Model = "HKY"
Params = 2.0, 1.0 
Freqs = 0.3, 0.2, 0.2, 0.3

## Root Section ################################################################
#
# Simulate a sequence that is 1000 nt long.

[Root]
Length = 1000

## Sim Section #################################################################
#
# Simulate 10 alignments

[Sim]
Reps = 10
"""

puppy = pd.PyDawg(simulationSeed=1111111, 
    inputFile=trickString, outputFile='basic_dna_cython_1111111.fasta', 
    [Segment(
        name='_default_',
        inherits_from='_initial_',
        output_rna=False,
        output_lowercase=False,
        output_keepempty=False,
        output_markins=True,
        substitution_model='jc',
        substitution_parameters=[0.3, 0.2],
        substitution_frequencies=[],
        substitution_rate_model='',
        substitution_rate_parameters=[],
        indel_model_insertion='POWER-LAW',
        indel_parameters_insertion=[0.3, 10.0],
        indel_rate_insertion='',
        indel_max_insertion='',
        indel_model_deletion='',
        indel_parameters_deletion='',
        indel_rate_deletion='',
        indel_max_deletion='',
        tree_model='',
        tree_parameters='',
        tree_tree='((A:0.1, B:0.4)C);', # Newick format, I assume
        tree_scale='',
        root_length=0,
        root_sequence='',
        root_rates=[],
        root_code=0,
        root_segment=0,
        root_gapoverlap=False)]
    )
puppy.run() # call walk from dawgcpp
puppy.printAlignments()

# segmentMap = {
#     "segment1" : pd.Segment(inheritsFrom="__default__",
# 	treeModel="((A:0.1, B:0.4)C);",
# 	substitutionModel="jc",
# 	substitutionParameters=[0.3, 0.2],
# 	indelModelInsertion="POWER-LAW",
# 	indelParametersInsertion= [0.3, 10.0],
# 	indelRateInsertion=0.335,
# 	indelMaxInsertion=10.0,
# 	indelModelDeletion="ZERO",
# 	indelParametersDeletion= [0.01, 4.0],
# 	indelRateDeletion=0.03322,
# 	indelMaxDeletion=4.0,
# 	rootCode=1,
# 	rootSequence="AATTTGGGGAAAAAAAAATTCC"),
#     "segment2" : pd.Segment(
# 	inheritsFrom="segment1",
# 	treeModel="((A:0.4, B:0.6)C);",
# 	substitutionModel="GAMMA",
# 	substitutionParameters=[0.3, 0.2],
# 	indelModelInsertion="POWER-LAW",
# 	indelParametersInsertion= [0.3, 10.0],
# 	indelRateInsertion=0.335,
# 	indelMaxInsertion=10.0,
# 	indelModelDeletion="ZERO",
# 	indelParametersDeletion= [0.01, 4.0],
# 	indelRateDeletion=0.03322,
# 	indelMaxDeletion=4.0,
#     rootCode=1,
# 	rootSequence="AATTTGGGGAAAAAAAAATTCC")}

    # # alignmentList is a list of strings containing the simulated alignments
    # alignmentList = dawg.fetchAlignments()

    # for alignment in alignmentList:
    #     seq = ''
    #     if builder.outputRna is False:
    #         seq = (alignment, DNAAlphabet())
    #     else:
    #         seq = (alignment, RNAAlphabet())

    #     seq.id = alignment.label
    #     seqRecord = SeqIO.parse(seq, builder.outputFormat)
    #     print("%s %s..." % (seqRecord.id, seqRecord.description[:50]))
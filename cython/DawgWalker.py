#!/usr/bin/python3

from os import sys
import PyDawg as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

# if __name__ == '__main__':
#     print('Welcome to the DawgBuilder demo')

trickString = \
"""
# Example: Specifying a recombination event and using autonamed sections
[[-]]
Subst.Model = jc
Root.Segment = 1
Root.Length = 60
Tree.Tree = ((A:0.02,B:0.02):0.2,(C:0.02):0.2);

[[-]]
Root.Code = 1
Root.Segment = 0
Root.Length = 60
Tree.Tree = ((A:0.02):0.2,(B:0.02,C:0.02):0.2);

[[-]]
Root.Segment = 2
"""

trickString2 = \
"""
# Example: Simulate DNA evolution along a tree
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

puppy = pd.PyDawg(output_file='basic_dna_cython_1111111.fasta', 
    simulation_seed=1111111, simulation_reps=10)
puppy.addSegment(
    name='_default_',
    inherits_from='_initial_',
    substitution_model='jc',
    substitution_parameters='0.3, 0.2',
    substitution_frequencies='',
    substitution_rate_model='',
    substitution_rate_parameters='',

    indel_model_insertion='POWER-LAW',
    indel_parameters_insertion='0.3, 10.0',
    indel_rate_insertion='',
    indel_max_insertion=1,
    indel_model_deletion='',
    indel_parameters_deletion='',
    indel_rate_deletion='',
    indel_max_deletion=1,

    tree_model='',
    tree_parameters='',
    tree_tree='\"((A:0.1, B:0.4)C);\"',
    tree_scale=0,

    root_length=0,
    root_sequence='',
    root_rates='',
    root_code=0,
    root_segment=0,
    root_gapoverlap=False,

    output_rna=False,
    output_lowercase=False,
    output_keepempty=False,
    output_markins=True)
puppy.addSegment(
    name='segment2',
    inherits_from='segment1',
    substitution_model='jc',
    substitution_parameters='0.3, 0.2',
    substitution_frequencies='',
    substitution_rate_model='',
    substitution_rate_parameters='',

    indel_model_insertion='POWER-LAW',
    indel_parameters_insertion='0.25, 10.0',
    indel_rate_insertion='',
    indel_max_insertion='',
    indel_model_deletion='',
    indel_parameters_deletion='',
    indel_rate_deletion='',
    indel_max_deletion='',

    tree_model='',
    tree_parameters='',
    tree_tree='((X:0.1, Y:0.4)Z);',
    tree_scale=0,

    root_length=0,
    root_sequence='',
    root_rates='',
    root_code=1,
    root_segment=0,
    root_gapoverlap=False,

    output_rna=False,
    output_lowercase=False,
    output_keepempty=False,
    output_markins=True)

puppy.echoSegments()
# puppy.run() # call walk from dawgcpp
# puppy.printAlignments()

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
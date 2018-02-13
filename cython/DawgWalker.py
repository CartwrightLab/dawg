#!/usr/bin/python3

from os import sys
import PyDawg as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

# First approach to doing this, we can just pass in a string to dawg,
# and Dawg Trick can parse just like a file
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
# Specify global options and the in/out
# puppy = pd.PyDawg(input=trickString, output='cython_trickstring_unaligned.fasta',
    # simulation_seed=1111111, simulation_reps=10)

# Better method that elimates parsing on Dawg's end,
# we configure the model arguments directly
goldenRetriever = pd.PyDawg(output='cython_modelargs.fasta',
    simulation_seed=1111111, simulation_reps=42)
goldenRetriever.addModelArgument(name='whatever', inherits_from='LUCA',
    tree="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
    substitution_model='HKY',
    substitution_params='2.0, 1.0',
    substitution_freqs='0.3, 0.2, 0.2, 0.3',
    root_length=1000)
goldenRetriever.bark()
# And then just run (configure, walk, aln)
# puppy.run() # call walk from dawgcpp
# puppy.printAlignments()

# We can also try elimating DAWG's output printer
# and use BioPython multiple-sequence alignment and fasta formatter

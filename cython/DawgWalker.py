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
# (Note to self, inheritance isn't working at the moment, and
# default section requires a tree?)
goldenRetriever = pd.PyDawg(output='pydawg_segments_1111111_10.fasta',
    simulation_seed=1111111, simulation_reps=10)
goldenRetriever.addModelArgument(name='stuff',
    tree="(M:0.80, N:0.20)P;",
    substitution_model="gtr",
    substitution_freqs='0.2, 0.3, 0.3, 0.2',
    substitution_params='2.0, 1.0, 3.0, 1.0, 1.0, 1.0',
    substitution_rate_model='zero',
    output_rna=True,
    output_markins=False,
    root_segment=0,
    root_length=10)
goldenRetriever.addModelArgument(name='whatever',
    tree="(A:0.2, B:0.3)C;",
    substitution_model='JC',
    root_length=10,
    root_segment=1)
goldenRetriever.addModelArgument(name='Coding',
    tree="(D:0.6, E:0.3)~F;",
    substitution_model='codmg',
    substitution_params='0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0',
    substitution_freqs='0.2, 0.3, 0.3, 0.2',
    root_length=20,
    root_segment=3)
goldenRetriever.addModelArgument(name='Noncoding2',
    tree='(x:0.3, y:0.1)C;',
    indel_params_insertion='1.03900, 100',
    indel_model_insertion="GEO",
    indel_rate_insertion='0.1',
    indel_max_insertion=100,
    indel_params_deletion='2, 10',
    indel_model_deletion="POWER-LAW",
    indel_rate_deletion='0.01',
    indel_max_deletion=10,
    root_length=10,
    root_sequence='acccctggggttaccccccccc',
    root_segment=4,
    output_rna=True)
# goldenRetriever.bark()
# And then just run (configure, walk, aln)
goldenRetriever.configureMatic()
goldenRetriever.walk()
goldenRetriever.write()

# We can also try elimating DAWG's output printer
# and use BioPython multiple-sequence alignment and fasta formatter

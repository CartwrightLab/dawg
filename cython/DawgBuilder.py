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

puppy = pd.PyDawg(simulationSeed=1111, simulationReps=222, 
    inputFile=trickString, outputFile='fasta:-', 
    outputRna=False, outputLowercase=False, outputKeepEmpty=False, outputMarkins=True)
puppy.run()


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
# for key, value in segmentMap.items():
#     print(key, 'corresponds to: ', segmentMap[key])
#     for segment in segmentMap[key]:
#         print()
# puppy.run()

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
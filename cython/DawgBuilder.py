#!/usr/bin/python3

from os import sys
import PyDawg as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

if __name__ == '__main__':
    print('Welcome to the DawgBuilder demo')

    afghaniHound = pd.PyDawg(
        simulationSeed=1111,
        simulationReps=222,
        inputFile='', outputFile='fasta:-',
        outputRna=False, outputLowercase=False, outputKeepEmpty=False,
        outputMarkins=True,
    	map(
    		("segment1", Segment(
        		inheritsFrom="__default__",
        		treeModel="((A:0.1, B:0.4)C);",
        		substModel="jc",
        		substParams=[0.3, 0.2],
        		indelInsModel="POWER-LAW",
        		indelInsParams= [0.3, 10.0],
        		indelInsRate=0.335,
        		indelInsMax=10.0,
        		indelDelModel="ZERO",
        		indelDelParams= [0.01, 4.0],
        		indelDelRate=0.03322,
        		indelDelMax=4.0
        		rootCode=1,
        		rootSequence="AATTTGGGGAAAAAAAAATTCC",
    		)),
            ("segment2",
    		inheritsFrom="segment1",
    		treeModel="((A:0.4, B:0.6)C);",
    		substModel="GAMMA",
    		substParams=[0.3, 0.2],
    		indelInsModel="POWER-LAW",
    		indelInsParams= [0.3, 10.0],
    		indelInsRate=0.335,
    		indelInsMax=10.0,
    		indelDelModel="ZERO",
    		indelDelParams= [0.01, 4.0],
    		indelDelRate=0.03322,
    		indelDelMax=4.0
    		rootCode=1,
    		rootSequence="AATTTGGGGAAAAAAAAATTCC")))

    # dawg = pd.PyDawg(dawgBuilder)
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
    

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

puppy = pd.PyDawg(input_file=trickString, output_file='cython_trickstring_unaligned.fasta',
    simulation_seed=1111111, simulation_reps=10)

# puppy.echoSegments()
puppy.run() # call walk from dawgcpp
# puppy.printAlignments()

# get unaligned sequences from DAWG
# Need to come up with a way of encoding the name of the evolved sequence here
evolvedSequenceList = [seq.split(":") for seq in puppy.fetchEvolvedSequences()]
for seqId, sequence in evolvedSequenceList.iteritems():
    # Assume DNA sequences for now
    seq = (sequence, DNAAlphabet())
    seq.id = seqId
    seqRecord = SeqIO.parse(seq, 'fasta')
    print("%s %s..." % (seqRecord.id, seqRecord.description[:50]))

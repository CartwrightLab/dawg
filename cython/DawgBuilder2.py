import PyDawg
import Bio
from os import sys

# PyDawg is a wrapper for DAWG and BioPython,
# accelerated by Cython. It interfaces with Dawg's C++ API where
# discrete sampling is done, and also wraps BioPython functionality

if __name__ == '__main__':
    # When a dog object is created, input can be a file or string
    # and output is simply a label which can end up being a file
    # or a string which can be printed to stdout
    # If no parameters are specified, then the only really functionality
    # that DAWG has is a RNG, or it makes assumptions about what letters
    # to use for the sequence alphabet with default parameters
    myDog = PyDawg.PyDawg()
    myDog.seedRng(seed=332)
    itr = myDog.generateSequence(type='dna')
    s = next(itr) # Using generator expressions to generate sequences incrementally
    # myDog.printSequence(format='fasta')

    # There are three ways to specify the parameters
    # for generating a sequence:
    # 1) Using named parameters
    # 2) Using an input file (DAWG file, which is MD-ish)
    # 3) Using a map {'seq_type': 'dna', ... }
    s = myDog.generateSequence(seq_type='dna',
        format='fasta', reps=10, seed=434)

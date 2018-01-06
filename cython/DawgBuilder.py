#!/usr/bin/python3

from os import sys
import getopt, mistune

import PyDawg as pd

# dawgArguments = [argument for string in sys.argv]
# print(dawgArguments)

renderer = mistune.Renderer(escape=True, hard_wrap=True)
markdown = mistune.Markdown(renderer=renderer)

basicDnaStr = """
# Tree
 - Tree = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);"

# Subst
 - Model = "HKY"
 - Params = 2.0, 1.0
 - Freqs = 0.3, 0.2, 0.2, 0.3

# Root
 - Length = 1000

# Sim
 - Reps = 10
"""

basicDnaOutStr = "fasta"

def parseDawgStr(dawgStr):
    parsed = markdown(dawgStr)
    print(parsed)

"""
example: ./DawgBuilder.py -s 444 -r 10 -i basicdna.dawg -o basicdna.fasta
"""
def parseDawgFile(dawgInput):
    dawgStr = open(dawgInput, "r")
    # parsedDawgStr = markdown(dawgStr)
    print(dawgStr.readlines())

def prepareDawgOutput(dawgOutput):
    pass

def main(argv):
    inputFile = ''
    outputFile = ''
    seed = ''
    reps = ''
    try:
        opts, args = getopt.getopt(argv, "hdr:s:i:o:", ["--help", "--debug", "--seed", "--reps", "--input","--output"])
    except getopt.GetoptError:
        print('error in arguments')
        # pd.PyDawg().help()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('help message')
            # pd.PyDawg().help()
            sys.exit()
        elif opt in ('-d', '--debug'):
            print('debug message')
            sys.exit()
        elif opt in ("-i", "--input"):
            inputFile = arg
        elif opt in ("-o", "--output"):
            outputFile = arg
        elif opt in ("-r", "--reps"):
            reps = arg
        elif opt in ("-s", "--seed"):
            seed = arg

    # if inputFile is '':
    #     pd.PyDawg().help()
    #     sys.exit(2)

    # print('Input file is: ', inputFile)
    # print('Output file is: ', outputFile)
    # print('reps count: ', reps)
    # print('seed count: ', seed)

    parseDawgFile(inputFile)

if __name__ == '__main__':
    print("Hello DAWG, I'm a Python")
    # for arg in sys.argv:
        # print(arg)
    # main(sys.argv[1:])
    parseDawgStr(basicDnaStr)

# The original example I got working
# Must use bytes on the string from Cython -> CPP (Cython book)
# donovan = pd.PyDawg(b"../examples/basic-dna.dawg", b"fasta:-", 10, 212121)
# donovan.run()

# class Bone:
#     __init__(self, bone):
#         self.bone = bone
#
# class Walk:
#     __init__(self, walk):
#         self.walk = walk
#
# class Segments:
#     __init__(self, name, inheritsFrom, *args):
#         self.header   = header
#         self.segments = args
#
# class Trick:
#     __init__(self, name, segments):
#         self.name     = name
#         self.segments = segments
#
# class DawgBuilder:
#     __init__(self, trick, walk, bone):
#         self.trick = trick
#         self.walk  = walk
#         self.bone  = bone

# PyDawg constructor can take in Tricks, Walks, and Bones
# akita = pd.pd(
# DawgBuilder(
#     Trick(name="BasicDnaExample",
#         Segments(name="__default__", inheritsFrom="LUCA",
#             (tree_tree="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
#             root_length=1000, subst_model='hky',
#             subst_freqs=[0.2, 0.3, 0.3, 0.2],
#             subst_params=[0.2, 1.0], sim_reps=10))),
#     Walk(), Bone("fasta:-")))

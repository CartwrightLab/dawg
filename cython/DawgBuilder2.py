import PyDawg
import Bio
from os import sys
import getopt, mistune

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

# # dawgArguments = [argument for string in sys.argv]
# # print(dawgArguments)
#
# renderer = mistune.Renderer(escape=True, hard_wrap=True)
# markdown = mistune.Markdown(renderer=renderer)
#
# basicDnaStr = """
# # Tree
#  - Tree = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);"
#
# # Subst
#  - Model = "HKY"
#  - Params = 2.0, 1.0
#  - Freqs = 0.3, 0.2, 0.2, 0.3
#
# # Root
#  - Length = 1000
#
# # Sim
#  - Reps = 10
# """
#
# basicDnaOutStr = "fasta"
#
# def parseDawgStr(dawgStr):
#     parsed = markdown(dawgStr)
#     print(parsed)
#
# """
# example: ./DawgBuilder.py -s 444 -r 10 -i basicdna.dawg -o basicdna.fasta
# """
# def parseDawgFile(dawgInput):
#     dawgStr = open(dawgInput, "r")
#     # parsedDawgStr = markdown(dawgStr)
#     print(dawgStr.readlines())
#
# def prepareDawgOutput(dawgOutput):
#     pass
#
# def main(argv):
#     inputFile = ''
#     outputFile = ''
#     seed = ''
#     reps = ''
#     try:
#         opts, args = getopt.getopt(argv, "hdr:s:i:o:", ["--help", "--debug", "--seed", "--reps", "--input","--output"])
#     except getopt.GetoptError:
#         print('error in arguments')
#         # pd.PyDawg().help()
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt in ('-h', '--help'):
#             print('help message')
#             # pd.PyDawg().help()
#             sys.exit()
#         elif opt in ('-d', '--debug'):
#             print('debug message')
#             sys.exit()
#         elif opt in ("-i", "--input"):
#             inputFile = arg
#         elif opt in ("-o", "--output"):
#             outputFile = arg
#         elif opt in ("-r", "--reps"):
#             reps = arg
#         elif opt in ("-s", "--seed"):
#             seed = arg
#
#     # if inputFile is '':
#     #     pd.PyDawg().help()
#     #     sys.exit(2)
#
#     # print('Input file is: ', inputFile)
#     # print('Output file is: ', outputFile)
#     # print('reps count: ', reps)
#     # print('seed count: ', seed)
#
#     parseDawgFile(inputFile)

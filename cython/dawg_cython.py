#!/usr/bin/python3
from os import sys
import PyDawg
# import pandas as pd

class TrickSection:
    pass

class DawgParams:
    pass

# Section commands are redundant
# Global commands are the same for all trick
# trickCommands = [
#     # Section 1
#     TrickSection(subst.model="",
#         subst.params="",
#         subst.freqs="",
#         subst.rate="",),
#     TrickSection(),
#     TrickSection(),
# ]

def ff():
    print("inside ff()")

if __name__ == '__main__':
    dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:/dev/null", 10, 212121)
    # dog1.run()
    # print("Stats for basic-dna.dawg trick file")
    # dog1.trickStats()
    #
    # numbers = []
    # for i in range(1000):
    #     numbers.append(dog1.rand(0, 100))
    # s = pd.Series(numbers)
    # s.value_counts()
    # df = pd.DataFrame(data=numbers)
    # print(df)
    # pd.DataFrame.hist(data=numbers)

    ts = TrickSection()
    ts.name = 'nothing'
    ts.val  = 'something'
    ts.thing = 'thing'
    ts.func = ff
    print(ts.func())

# [REGULAR PARAMETERS]
# Subst.Model - The identifier of the substitution model, e.g. JC, GTR, WAG,
#   CODGY.
# Subst.Params - A list specifying the parameters of the substitution model.
#   Model Dependant.
# Subst.Freqs - A list specifying the stationary frequencies of nucleotides,
#   amino acids, or codons. Model Dependant.
# Subst.Rate.Model - The identifier of the heterogeneous rate model, e.g.
#   CONST, GAMMA, or ZERO.
# Subst.Rate.Params - The parameters of the rate model.  Model Dependant.
# Indel.Model.Ins - The identifiers of the insertion models, e.g. USER, GEO,
#   POWER-LAW.
# Indel.Params.Ins - The parameters of the insertion models.  Model Dependant.
# Indel.Rate.Ins - The per-substitution rates of the mixture of insertion models.
# Indel.Max.Ins - The maximum size of an insertion
# Indel.Model.Del - The identifiers of the deletion models, e.g. USER, GEO,
#   POWER-LAW.
# Indel.Params.Del - The parameters of the deletion models.  Model Dependant.
# Indel.Rate.Del - The per-substitution rates of the mixture of deletion models.
# Indel.Max.Del - The maximum size of a deletion.
# Tree.Model - The identifier of the tree model.
# Tree.Params - The parameters of the tree model.  Model Dependant.
# Tree.Tree - The tree or tree template.
# Tree.Scale - Branch-lengths are scaled by this number in the simulation.
# Root.Length - The length of a randomly generated root sequence.
# Root.Seq - A specific root sequence.
# Root.Rates - The heterogeneous rates of the root sequence.
# Root.Code - The genetic code used when simulating codon evolution.
# Root.Segment - The segment number that the root belongs too.
# Root.Gapoverlap - Allow upstream deletions to affect this segment.
# Output.Markins - Distinguish insertions from deletions.
# Output.Keepempty - Keep empty columns instead of deleting them in the alignment.
# Output.Lowercase - Use lowercase for sequence output.
# Output.Rna - Output an RNA sequence instead of a DNA sequence
#
# [GLOBAL PARAMETERS]
# Output.Block.Head - Text that will be written to the beginning of output.
# Output.Block.Tail - Text that will be written to the end of output.
# Output.Block.Before - Text that will be written before every replicate.
# Output.Block.After - Text that will be written after every replicate.
# Output.Block.Between - Text that will be written between replicates.
# Output.File - Path to the output file.
# Output.Split - Output each replicate to its own file.
# Output.Append - Append results to existing file.
# Output.Label - label each simulation with a unique id.
# Sim.Reps - Number of simulation replicates.
# Sim.Seed - The seed of the random number generator

    params = DawgParams()
    params.SubstitutionModel            = 'jc'
    params.SubstitutionParameters       = 'aaa'
    params.SubstitutionFrequencies      = 'aa'
    params.SubstitutionRateModel        = 'bb'
    params.SubstitutionRateParameters   = 'cc'
    params.IndelParametersInsertion     = 'dd'
    params.IndelRateInsertion           = '1'
    params.IndelMaximumInsertion        = '2'
    params.IndelMinimumInsertion        = '22'
    params.IndelModelInsertion          = '3'
    params.IndelModelDeletion           = '33'
    params.IndelParametersDeletion      = '4'
    params.IndelRateDeletion            = '5'
    params.IndelMinimumDeletion         = '55'
    params.IndelMaximumDeletion         = '6'
    params.IndelDeletionCharacter       = 'd'
    params.IndelInsertionCharacter      = 'i'
    params.TreeModel                    = 'juniper'
    params.TreeParameters               = '3343'
    params.Tree                         = '((sajsjajaj): 0303);'
    params.TreeScale                    = '100%'
    params.RootLength                   = 1000
    params.RootSequence                 = 'ADDDTADDGGCCCT'
    params.RootRates                    = 300
    params.RootCode                     = 100
    params.RootSegment                  = 1
    params.RootGapOverlap               = 't'
    params.SimulationRepititions        = 3
    params.SimulationSeed               = 4444
    params.SimulationDawgFile           = 'akita.dawg'
    params.SimulationFormat             = 'fasta'
    params.SimulationOutputFile         = 'akita.fasta'

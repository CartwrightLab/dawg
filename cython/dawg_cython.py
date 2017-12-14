#!/usr/bin/python3
from os import sys
# import PyDawg

# dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:/dev/null", 10, 212121)
# dog1.run()
# print("Stats for basic-dna.dawg trick file")
# dog1.trickStats()
# dog1.rand(0, 100))


class TrickSection:

    populateSection(self, name, inheritsFrom, map db):
        self.name = name
        self.inheritsFrom = inheritsFrom
        self.db = db

    compileNewickTree():
        pass

# The list of sections to breed
# There is a priority queue that sorts the sections by
# how easy they will be to evolve?
class DawgBreed:

    insert(self, section)
        self.pqueue.insert(section)

def basicdna():
    # setup the basic-dna example from dawg without parsing
    ts1 = TrickSection()
    ts1.populateSection(name="basicdna", inheritsFrom="LUCA", [
        "Tree.Tree"="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
        "Subst.Model"="HKY",
        "Subst.Params"="2.0, 1.0",
        "Subst.Freqs"="0.3, 0.2, 0.2, 0.3",
        "Root.Length"="1000",
        "Sim.Reps"="10"
    ])

    ts1.compileNewickTree()

def help():
    print("""
dawg 2-current-rUnknown
    Copyright (C) 2004-2013  Reed A. Cartwright, PhD <cartwright@asu.edu>

Usage:
  dawg [options] trick-1.dawg trick-2.dawg ...

Allowed Options:
  --version                    display version information
  --help-trick                 display description of common control variables
  --help                       display help message
  -o [ --output ] arg          output to this file
  --seed arg (=0)              PRNG seed
  --reps arg (=0)              the number of alignments to generate
  --split [=arg(=on)] (=null)  split output into separate files
  --append [=arg(=on)] (=null) append output to file
  --label [=arg(=on)] (=null)  label each simulation with a unique id
  --arg-file arg               read arguments from file

[REGULAR PARAMETERS]
Subst.Model - The identifier of the substitution model, e.g. JC, GTR, WAG,
  CODGY.
Subst.Params - A list specifying the parameters of the substitution model.
  Model Dependant.
Subst.Freqs - A list specifying the stationary frequencies of nucleotides,
  amino acids, or codons. Model Dependant.
Subst.Rate.Model - The identifier of the heterogeneous rate model, e.g.
  CONST, GAMMA, or ZERO.
Subst.Rate.Params - The parameters of the rate model.  Model Dependant.
Indel.Model.Ins - The identifiers of the insertion models, e.g. USER, GEO,
  POWER-LAW.
Indel.Params.Ins - The parameters of the insertion models.  Model Dependant.
Indel.Rate.Ins - The per-substitution rates of the mixture of insertion models.
Indel.Max.Ins - The maximum size of an insertion
Indel.Model.Del - The identifiers of the deletion models, e.g. USER, GEO,
  POWER-LAW.
Indel.Params.Del - The parameters of the deletion models.  Model Dependant.
Indel.Rate.Del - The per-substitution rates of the mixture of deletion models.
Indel.Max.Del - The maximum size of a deletion.
Tree.Model - The identifier of the tree model.
Tree.Params - The parameters of the tree model.  Model Dependant.
Tree.Tree - The tree or tree template.
Tree.Scale - Branch-lengths are scaled by this number in the simulation.
Root.Length - The length of a randomly generated root sequence.
Root.Seq - A specific root sequence.
Root.Rates - The heterogeneous rates of the root sequence.
Root.Code - The genetic code used when simulating codon evolution.
Root.Segment - The segment number that the root belongs too.
Root.Gapoverlap - Allow upstream deletions to affect this segment.
Output.Markins - Distinguish insertions from deletions.
Output.Keepempty - Keep empty columns instead of deleting them in the alignment.
Output.Lowercase - Use lowercase for sequence output.
Output.Rna - Output an RNA sequence instead of a DNA sequence

[GLOBAL PARAMETERS]
Output.Block.Head - Text that will be written to the beginning of output.
Output.Block.Tail - Text that will be written to the end of output.
Output.Block.Before - Text that will be written before every replicate.
Output.Block.After - Text that will be written after every replicate.
Output.Block.Between - Text that will be written between replicates.
Output.File - Path to the output file.
Output.Split - Output each replicate to its own file.
Output.Append - Append results to existing file.
Output.Label - label each simulation with a unique id.
Sim.Reps - Number of simulation replicates.
Sim.Seed - The seed of the random number generator
    """)

help()

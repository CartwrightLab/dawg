# cython: language_level = 3
# distutils: language = c++
# distutils: sources = dawg.cpp

from Bio import SeqIO

#from cpython.array cimport array
#from libcpp.list cimport list
#from libcpp.array cimport array
from libcpp.string cimport string
#from libcpp.map cimport map

cdef extern from "dawg.hpp" namespace "dawg":
    cdef cppclass Dawg:
        Dawg()
        Dawg(unsigned int)
        Dawg(string, string, unsigned int, unsigned int)
        void run()
        void bark()
        unsigned int rand(unsigned int, unsigned int)
        void trickStats()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, *args):
        if len(args) == 4:
            self._thisptr = new Dawg(args[0], args[1], args[2], args[3])
        elif len(args) == 1:
            self._thisptr = new Dawg(args[0])
        else:
            self._thisptr = new Dawg()

        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

    # Get the alignment string
    cpdef void bark(self):
        # self._thisptr.bark()

    cpdef unsigned int rand(self, a, b):
        return self._thisptr.rand(a, b)

    cpdef void trickStats(self):
        self._thisptr.trickStats()

    # IO functions
    # trick takes in a bunch of
    def trick(self, s, **kwargs):
        self.trickQueue = kwargs
        self.trickQueue.label = s

    def outputHandler(self):
        pass

    def help(self):
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

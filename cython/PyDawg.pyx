# cython: language_level = 3
# distutils: language = c++
# distutils: sources = dawg.cpp

#from cpython.array cimport array
#from libcpp.list cimport list
#from libcpp.array cimport array
from libcpp.string cimport string
#from libcpp.map cimport map

cdef extern from "dawg.hpp" namespace "dawg":
    cdef cppclass Dawg:
        # Dawg()
        # Dawg(unsigned int)
        # Dawg(string, string, unsigned int, unsigned int)

        # del model, del rate, del param 1, del param 2, max del
        Dawg(string, float, float, float, float)
        void run()
        void bark()
        unsigned int rand(unsigned int, unsigned int)
        void trickStats()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self,
        simulationSeed=0,
        simulationReps=10,
        inputFile='',
        outputFile='fasta:-',
        outputSplit=False,
        outputAppend=False,
        outputLabel=False,
        outputRna=False,
        outputProtein=False,
        outputLowercase=False,
        outputKeepEmpty=False,
        outputMarkins=False,
        segmentMap={}):

        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

    # Security?
    cpdef void bark(self):
        pass

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
        print("Mind your own business")

class Segment:

    def __init__(self,
        name='',
        inheritsFrom='',
        substitutionModel='',
        substitutionParameters=[],
        substitutionFrequencies=[],
        substitutionRateModel='',
        substitutionRateParameters=[],
        indelModelInsertion='',
        indelParametersInsertion='',
        indelRateInsertion='',
        indelMaxInsertion='',
        indelModelDeletion='',
        indelParametersDeletion='',
        indelRateDeletion='',
        indelMaxDeletion='',
        treeModel='',
        treeParameters='',
        treeTree='', # Newick format, I assume
        treeScale='',
        rootLength=0,
        rootSequence='',
        rootRates=[],
        rootCode=0,
        rootSegment=0,
        rootGapOverlap=False):
        
        self.name = name
        self.inheritsFrom = inheritsFrom
        
        self.substitutionModel = substitutionModel
        self.substitutionParameters = substitutionParameters
        self.substitutionFrequencies = substitutionFrequencies
        self.substitutionRateModel = substitutionRateModel
        self.substitutionRateParameters = substitutionRateParameters

        self.indelModelInsertion = indelModelInsertion
        self.indelParametersInsertion = indelParametersInsertion
        self.indelRateInsertion = indelRateInsertion
        self.indelMaxInsertion = indelMaxInsertion
        
        self.indelModelDeletion = indelModelDeletion
        self.indelParametersDeletion = indelParametersDeletion
        self.indelRateDeletion = indelRateDeletion
        self.indelMaxDeletion = indelMaxDeletion

        self.treeModel = treeModel
        self.treeParameters = treeParameters
        self.treeTree = treeTree
        self.treeScale = treeScale

        self.rootLength = rootLength
        self.rootSequence = rootSequence
        self.rootRates = rootRates
        self.rootCode = rootCode
        self.rootSegment = rootSegment
        self.rootGapOverlap = rootGapOverlap
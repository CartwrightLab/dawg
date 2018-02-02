# cython: language_level = 3
# distutils: language = c++
# distutils: sources = dawg.cpp

#from cpython.array cimport array
#from libcpp.list cimport list
#from libcpp.array cimport array
from libcpp.string cimport string
# from libcpp.map cimport map as cppMap
# from libcpp.vector cimport vector

"""
This is Cython's way of interfacing with C++
"""
cdef extern from "dawg.hpp" namespace "dawg":
    cdef cppclass Dawg:
        Dawg()
        Dawg(unsigned int seed)
        # Dawg(cppMap, string, unsigned int, unsigned int)
        # infile, outfile, seed, reps
        Dawg(string, string, unsigned int, unsigned int)

        void run()
        void printAlignments()
        void bark()
        unsigned int rand(unsigned int, unsigned int)
        void trickStats()

"""
The PyDawg module is the official Python module we can import,
it acts as an opaque pointer to the Cython cppclass
"""
cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self,
        simulation_seed=0, # global option
        simulation_reps=10, # global option
        input_file='', # global option
        output_file='fasta:-', # global option
        output_split=False, # global option
        output_append=False, # global option
        output_label=False, # global option
        # model options
        segments=[]):

        self.simulation_seed = simulation_seed
        self.simulation_reps = simulation_reps
        self.input_file = input_file
        self.output_file = output_file
        self.output_split = output_split
        self.output_append = output_append
        self.output_label = output_label
    
        self.segments = segments

        # Input file is the main variant which we use 
        # to determine CPP Dawg args
        # If the input is empty, then look at segments
        if not inputFile:

        else:
            self._thisptr = new Dawg(
                inputFile.encode(), outputFile.encode(), 
                simulationSeed, simulationReps)

        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

    cpdef void printAlignments(self):
        self._thisptr.printAlignments()

    # Security?
    cpdef void bark(self):
        pass

    cpdef unsigned int rand(self, a, b):
        return self._thisptr.rand(a, b)

    def help(self):
        print("Mind your own business")

class Segment:

    def __init__(self,
        name='_default_',
        inherits_from='_initial_',
        output_rna=False,
        # outputProtein=False,
        output_lowercase=False,
        output_keepempty=False,
        output_markins=False,
        substitution_model='',
        substitution_parameters=[],
        substitution_frequencies=[],
        substitution_rate_model='',
        substitution_rate_parameters=[],
        indel_model_insertion='',
        indel_parameters_insertion='',
        indel_rate_insertion='',
        indel_max_insertion='',
        indel_model_deletion='',
        indel_parameters_deletion='',
        indel_rate_deletion='',
        indel_max_deletion='',
        tree_model='',
        tree_parameters='',
        tree_tree='', # Newick format, I assume
        tree_scale='',
        root_length=0,
        root_sequence='',
        root_rates=[],
        root_code=0,
        root_segment=0,
        root_gapoverlap=False):
        
        self.name = name
        self.inherits_from = inherits_from
        
        self.output_rna = output_rna
        self.output_protein = output_protein
        self.output_lowercase = output_lowercase
        self.output_keepempty = output_keepempty
        self.output_markins = output_markins

        self.substitution_model = substitution_model
        self.substitution_parameters = substitution_parameters
        self.substitution_frequencies = substitution_frequencies
        self.substitution_rate_model = substitution_rate_model
        self.substitution_rate_parameters = substitution_rate_parameters

        self.indel_model_insertion = indel_model_insertion
        self.indel_parameters_insertion = indel_parameters_insertion
        self.indel_rate_insertion = indel_rate_insertion
        self.indel_max_insertion = indel_max_insertion
        
        self.indel_model_deletion = indel_model_deletion
        self.indel_parameters_deletion = indel_parameters_deletion
        self.indel_rate_deletion = indel_rate_deletion
        self.indel_max_deletion = indel_max_deletion

        self.tree_model = tree_model
        self.tree_parameters = tree_parameters
        self.tree_tree = tree_tree
        self.tree_scale = tree_scale

        self.root_length = root_length
        self.root_sequence = root_sequence
        self.root_rates = root_rates
        self.root_code = root_code
        self.root_segment = root_segment
        self.root_gapoverlap = root_gapoverlap
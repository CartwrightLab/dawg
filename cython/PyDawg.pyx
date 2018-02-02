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
        # outfile, seed, reps
        Dawg(string, unsigned int, unsigned int)
        # infile, outfile, seed, reps
        Dawg(string, string, unsigned int, unsigned int)

        void addSegment(string name,
            # string inherits_from,
            string substitution_model,
            string substitution_parameters,
            string substitution_frequencies,
            string substitution_rate_model,
            string substitution_rate_parameters,

            string indel_model_insertion,
            string indel_parameters_insertion,
            string indel_rate_insertion,
            unsigned int indel_max_insertion,
            string indel_model_deletion,
            string indel_parameters_deletion,
            string indel_rate_deletion,
            unsigned int indel_max_deletion,

            string tree_model,
            string tree_parameters,
            string tree_tree,
            double tree_scale,

            unsigned int root_length,
            string root_sequence,
            string root_rates,
            unsigned int root_code,
            unsigned int root_segment,
            bint root_gapoverlap,

            bint output_rna,
            bint output_lowercase,
            bint output_keepempty,
            bint output_markins)
        void echoSegments()
        void run()
        void printAlignments()
        void bark()
        unsigned int rand(unsigned int, unsigned int)
        void trickStats()

"""
The PyDawg module is the official Python module we can import,
it acts as an opaque pointer to the Cython cppclass
The constructor takes in the global options of the simulation
"""
cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, input_file='', output_file='fasta:-', 
        simulation_seed=0, simulation_reps=10):
        # output_split=False,
        # output_append=False,
        # output_label=False):
        # model options
        # segments=[]):

        if not input_file:
            self._thisptr = new Dawg(output_file.encode(), simulation_seed, simulation_reps)
        else:
            self._thisptr = new Dawg(input_file.encode(), output_file.encode(), 
                simulation_seed, simulation_reps)

        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    # encode all strings by default, utf-8
    cpdef void addSegment(self,
        name='_default_',
        inherits_from='_initial_',
        substitution_model='jc',
        substitution_parameters='0.3, 0.2',
        substitution_frequencies='',
        substitution_rate_model='',
        substitution_rate_parameters='',

        indel_model_insertion='POWER-LAW',
        indel_parameters_insertion='0.3, 10.0',
        indel_rate_insertion='',
        indel_max_insertion=1,
        indel_model_deletion='',
        indel_parameters_deletion='',
        indel_rate_deletion='',
        indel_max_deletion=1,

        tree_model='',
        tree_parameters='',
        tree_tree='\"((A:0.1, B:0.4)C);\"',
        tree_scale=0,

        root_length=0,
        root_sequence='',
        root_rates='',
        root_code=0,
        root_segment=0,
        root_gapoverlap=False,

        output_rna=False,
        output_lowercase=False,
        output_keepempty=False,
        output_markins=True):
            self._thisptr.addSegment(name.encode(),
            #inheritsFrom.encode(), 
            substitution_model.encode(),
            substitution_parameters.encode(),
            substitution_frequencies.encode(),
            substitution_rate_model.encode(),
            substitution_rate_parameters.encode(),
            
            indel_model_insertion.encode(),
            indel_parameters_insertion.encode(),
            indel_rate_insertion.encode(),
            indel_max_insertion,
            indel_model_deletion.encode(),
            indel_parameters_deletion.encode(),
            indel_rate_deletion.encode(),
            indel_max_deletion,

            tree_model.encode(),
            tree_parameters.encode(),
            tree_tree.encode(),
            tree_scale,

            root_length,
            root_sequence.encode(),
            root_rates.encode(),
            root_code,
            root_segment,
            root_gapoverlap,

            output_rna, output_lowercase, output_keepempty, output_markins)

    cpdef void echoSegments(self):
        self._thisptr.echoSegments()

    cpdef void run(self):
        self._thisptr.run()

    cpdef void printAlignments(self):
        self._thisptr.printAlignments()

    # Wolves howl, dogs bark ... well sometimes they howl too, but wolves don't bark?
    cpdef void bark(self):
        pass

    cpdef unsigned int rand(self, a, b):
        return self._thisptr.rand(a, b)

    def help(self):
        print("Help me.")

# class Segment:

#     def __init__(self,
#         name='_default_',
#         inherits_from='_initial_',
#         output_rna=False,
#         # outputProtein=False,
#         output_lowercase=False,
#         output_keepempty=False,
#         output_markins=False,
#         substitution_model='',
#         substitution_parameters=[],
#         substitution_frequencies=[],
#         substitution_rate_model='',
#         substitution_rate_parameters=[],
#         indel_model_insertion='',
#         indel_parameters_insertion='',
#         indel_rate_insertion='',
#         indel_max_insertion='',
#         indel_model_deletion='',
#         indel_parameters_deletion='',
#         indel_rate_deletion='',
#         indel_max_deletion='',
#         tree_model='',
#         tree_parameters='',
#         tree_tree='', # Newick format, I assume
#         tree_scale='',
#         root_length=0,
#         root_sequence='',
#         root_rates=[],
#         root_code=0,
#         root_segment=0,
#         root_gapoverlap=False):
        
#         self.name = name
#         self.inherits_from = inherits_from
        
#         self.output_rna = output_rna
#         self.output_protein = output_protein
#         self.output_lowercase = output_lowercase
#         self.output_keepempty = output_keepempty
#         self.output_markins = output_markins

#         self.substitution_model = substitution_model
#         self.substitution_parameters = substitution_parameters
#         self.substitution_frequencies = substitution_frequencies
#         self.substitution_rate_model = substitution_rate_model
#         self.substitution_rate_parameters = substitution_rate_parameters

#         self.indel_model_insertion = indel_model_insertion
#         self.indel_parameters_insertion = indel_parameters_insertion
#         self.indel_rate_insertion = indel_rate_insertion
#         self.indel_max_insertion = indel_max_insertion
        
#         self.indel_model_deletion = indel_model_deletion
#         self.indel_parameters_deletion = indel_parameters_deletion
#         self.indel_rate_deletion = indel_rate_deletion
#         self.indel_max_deletion = indel_max_deletion

#         self.tree_model = tree_model
#         self.tree_parameters = tree_parameters
#         self.tree_tree = tree_tree
#         self.tree_scale = tree_scale

#         self.root_length = root_length
#         self.root_sequence = root_sequence
#         self.root_rates = root_rates
#         self.root_code = root_code
#         self.root_segment = root_segment
#         self.root_gapoverlap = root_gapoverlap
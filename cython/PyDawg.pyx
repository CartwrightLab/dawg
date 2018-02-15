# cython: language_level = 3
# distutils: language = c++
# distutils: sources = dawg.cpp

from libcpp.string cimport string

# Try to keep params as generic as possible,
# Let the user or the CPP-side handle validation
cdef extern from "dawg.hpp" namespace "dawg":
    cdef cppclass Dawg:
        Dawg()
        Dawg(unsigned int seed)
        Dawg(string input, string output,
            unsigned int seed, unsigned int reps)

        void addModelArgument(string name,
            string inherits_from,
            string substitution_model,
            string substitution_params,
            string substitution_freqs,
            string substitution_rate_model,
            string substitution_rate_params,

            string indel_model_insertion,
            string indel_params_insertion,
            string indel_rate_insertion,
            unsigned int indel_max_insertion,
            string indel_model_deletion,
            string indel_params_deletion,
            string indel_rate_deletion,
            unsigned int indel_max_deletion,

            string tree,
            string tree_model,
            string tree_params,
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
        void configureMatic()
        void walk()
        void write()
        string getAlignments()
        unsigned int rand(unsigned int, unsigned int)
        void bark()

"""
The PyDawg module is the official Python module we can import,
it acts as an opaque pointer to the Cython cppclass
The constructor takes in the global options of the simulation
Validation and string splitting happens on CPP side
"""
cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, input='', output='fasta:-',
        simulation_seed=0, simulation_reps=10):

        self._thisptr = new Dawg(input.encode(), output.encode(),
            simulation_seed, simulation_reps)

        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    # encode all strings by default, utf-8
    # default params are taken from src/include/dawg/ma.h
    cpdef void addModelArgument(self,
        name='_initial_',
        inherits_from='_default_',
        substitution_model='jc',
        substitution_params='',
        substitution_freqs='',
        substitution_rate_model='const',
        substitution_rate_params='',

        indel_model_insertion='user',
        indel_params_insertion='1',
        indel_rate_insertion='',
        indel_max_insertion=100,
        indel_model_deletion='user',
        indel_params_deletion='1',
        indel_rate_deletion='',
        indel_max_deletion=100,

        tree='',
        tree_model='user',
        tree_params='',
        tree_scale=1,

        root_length=0,
        root_sequence='',
        root_rates='', # this isn't being used but kept it
        root_code=0,
        root_segment=1,
        root_gapoverlap=True,

        output_rna=False,
        output_lowercase=False,
        output_keepempty=True,
        output_markins=False):
            self._thisptr.addModelArgument(
            name.encode(),
            inherits_from.encode(),
            substitution_model.encode(),
            substitution_params.encode(),
            substitution_freqs.encode(),
            substitution_rate_model.encode(),
            substitution_rate_params.encode(),

            indel_model_insertion.encode(),
            indel_params_insertion.encode(),
            indel_rate_insertion.encode(),
            indel_max_insertion,
            indel_model_deletion.encode(),
            indel_params_deletion.encode(),
            indel_rate_deletion.encode(),
            indel_max_deletion,

            tree.encode(),
            tree_model.encode(),
            tree_params.encode(),
            tree_scale,

            root_length,
            root_sequence.encode(),
            root_rates.encode(),
            root_code,
            root_segment,
            root_gapoverlap,

            output_rna,
            output_lowercase,
            output_keepempty,
            output_markins)

    cpdef void configureMatic(self):
        self._thisptr.configureMatic()

    cpdef void walk(self):
        self._thisptr.walk()

    cpdef void write(self):
        self._thisptr.write()

    # The return alignments are encrypted
    # example: "label1:seq1;label2:seq2;labelN:seqN"
    # no semi colon at the end!
    cpdef string getAlignments(self):
        return self._thisptr.getAlignments()

    cpdef unsigned int rand(self, a, b):
        return self._thisptr.rand(a, b)

    cpdef void bark(self):
        self._thisptr.bark()

    def help(self):
        print("TODO, add help message.")

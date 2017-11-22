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

    cpdef void bark(self):
        self._thisptr.bark()

    cpdef unsigned int rand(self, a, b):
        return self._thisptr.rand(a, b)

    cpdef void trickStats(self):
        self._thisptr.trickStats()

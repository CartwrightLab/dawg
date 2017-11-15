# cython: language_level = 3
# distutils: language = c++
# distutils: sources = dawg.cpp

#from cpython.array cimport array
#from libcpp.list cimport list
#from libcpp.array cimport array
from libcpp.string cimport string

cdef extern from "dawg.hpp" namespace "dawg":
    cdef cppclass Dawg:
        Dawg()
        Dawg(string, string, unsigned int, unsigned int)
        void run()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, *args):
        self._thisptr = new Dawg(args[0], args[1], args[2], args[3])
        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

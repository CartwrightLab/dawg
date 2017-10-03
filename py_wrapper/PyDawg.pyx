# distutils: language = c++
# distutils: sources = dawg.cpp

from cpython.array cimport array

cdef extern from "dawg.hpp" namespace "dawg":
    unsigned int N
    cdef cppclass Dawg:
        Dawg(int argc, char* argv)
        void run()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, *args):

        self._thisptr = new Dawg(len(args), array("B", args).data.as_chars)
        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

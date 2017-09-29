# distutils: language = c++
# distutils: sources = dawg.cpp

cdef extern from "dawg.hpp" namespace "dawg":
    unsigned int N
    cdef cppclass Dawg:
        Dawg()
        void run()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self):
        self._thisptr = new Dawg()
        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        print("Hello run method")

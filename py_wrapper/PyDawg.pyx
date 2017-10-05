# distutils: language = c++
# distutils: sources = dawg.cpp

from cpython.array cimport array
from libcpp.list cimport list
#from libcpp.array cimport array
from libcpp.string cimport string

cdef extern from "dawg.hpp" namespace "dawg":
    unsigned int N
    cdef cppclass Dawg:
        Dawg()
        Dawg(list[string] args)
        void run()

cdef class PyDawg:

    cdef Dawg *_thisptr

    def __cinit__(self, args):
        cdef list[string] = args
        # for arg in args:
        #     print("{}".format(arg))
        # print("type of args:{}".format(args))
        # print("type of arr:{}".format(arr))
        #self._thisptr = new Dawg(array("B", args).data.as_chars)
        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        self._thisptr.run()

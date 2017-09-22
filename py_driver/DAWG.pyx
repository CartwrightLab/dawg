# distutils: language = c++
# distutils: sources = dawg.cpp

cdef extern from "dawg_app.h":
    cdef cppclass dawg_app:
        dawg_app(int argc, char* argv[])
        void run()

cdef class DAWG:

    cdef DAWG *_thisptr

    def __cinit__(self, int argc, char* argv[]):
        self._thisptr = new dawg_app(argc, argv)
        if self._thisptr == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    cpdef void run(self):
        return self._thisptr.run()

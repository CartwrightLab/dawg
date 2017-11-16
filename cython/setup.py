from distutils.core import setup, Extension
from Cython.Build import cythonize

# Extension name likes to match pyx name according to an error msg
ext = Extension("PyDawg",
                sources=["PyDawg.pyx", "dawg.cpp", "../src/lib/output.cpp",
                    "../src/lib/parse.cpp", "../src/lib/mutt.cpp", "../src/lib/models.cpp",
                    "../src/lib/ma.cpp", "../src/lib/matic.cpp"],
                include_dirs=["../src/include"],
                # libraries=["dawg2"],
                # library_dirs=["/usr/local/lib"],
                # runtime_library_dirs=[],
                # extra_compile_args=["-fPIC"],
                # extra_link_args=["-static"],
                # extra_objects=["ma.o"],
                #["boost_program_options"],
                language="c++")

setup(name="DAWG",
      ext_modules=cythonize(ext))

from distutils.core import setup, Extension
from Cython.Build import cythonize

# Extension name likes to match pyx name according to an error msg
ext = Extension("PyDawg",
                sources=["PyDawg.pyx", "dawg.cpp"],
                language="c++")

setup(name="DAWG",
      ext_modules=cythonize(ext))

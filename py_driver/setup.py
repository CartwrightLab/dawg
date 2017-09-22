from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("DAWG",
                sources=["DAWG.pyx", "../src/dawg.cpp"],
                language="c++")

setup(name="DAWG",
      ext_modules=cythonize(ext))

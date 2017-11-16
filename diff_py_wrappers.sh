#!/usr/bin/bash

python boost.python/dawg.py > boostPython.fasta
python cython/DawgDriver.py > cythonDawg.fasta

diff boost.python/boostPython.fasta cython/cythonDawg.fasta
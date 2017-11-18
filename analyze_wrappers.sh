#!/usr/bin/bash

# I used the basic-dna.dawg for all these (10 reps, seed=212121)
BOOST_FASTA_FILE="boostPython.fasta"
CYTHON_FASTA_FILE="cython.fasta"
SWIG_FASTA_FILE="swig.fasta"
DEV_FASTA_FILE="dev.fasta"

# cd boost.python && python dawg_boost.py > $BOOST_FASTA_FILE
# cd ../
# cd cython       && python dawg_cython.py > $CYTHON_FASTA_FILE
# cd ../
# cd swig && python dawg_swig.py > $SWIG_FASTA_FILE
# cd ../

diff $DEV_FASTA_FILE boost.python/$BOOST_FASTA_FILE
diff $DEV_FASTA_FILE cython/$CYTHON_FASTA_FILE
diff $DEV_FASTA_FILE swig/$SWIG_FASTA_FILE

# # Time Cython (approx 3m51.320s)
# # size of PyDawg.so file is about 2gb
# cd cython
# for (( i = 0; i < 10000; i++ ))
# do
# 	python dawg_cython.py > /dev/null
# done
# cd ../

# # Time BoostPython (approx 6m45.847s)
# # size of dawg.so file is about 702mb
# cd boost.python
# for (( i = 0; i < 10000; i++ ))
# do
# 	python dawg_boost.py > /dev/null
# done
# cd ../

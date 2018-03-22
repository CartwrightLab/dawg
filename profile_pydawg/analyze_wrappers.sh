#!/usr/bin/bash

# like wrappers do

# The **basic-dna.dawg** trick file with 10 reps, seed=212121 is used
BOOST_FASTA_FILE="boostPython.fasta"
CYTHON_FASTA_FILE="cython.fasta"
SWIG_FASTA_FILE="swig.fasta"
DEV_FASTA_FILE="dev.fasta"

# cd boost.python && python dawg_boost.py > $BOOST_FASTA_FILE
# cd ../
# cd cython       && python dawg_cython.py > $CYTHON_FASTA_FILE
# cd ../
# cd swig		  && python dawg_swig.py > $SWIG_FASTA_FILE
# cd ../
#
# diff $DEV_FASTA_FILE boost.python/$BOOST_FASTA_FILE
# diff $DEV_FASTA_FILE cython/$CYTHON_FASTA_FILE
# diff $DEV_FASTA_FILE swig/$SWIG_FASTA_FILE

# arg1=dir, arg2=filename
function profile {
	cd $1
	for (( i = 0; i < 10000; i++ ))
	do
		python $2 > /dev/null
	done
	cd ../
}

# Time SWIG (approx 7m17.167s)
# size of dawg.so file is about 6.6mb
# # Time BoostPython (approx 6m45.847s)
# # size of dawg.so file is about 702mb
# # Time Cython (approx 3m51.320s)
# # size of PyDawg.so file is about 2gb
time profile swig/ dawg_swig.py
time profile cython/ dawg_cython.py
time profile boost.python/ dawg_boost.py

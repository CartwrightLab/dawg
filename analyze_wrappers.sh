#!/usr/bin/bash

BOOST_FASTA_FILE="boostPython.fasta"
CYTHON_FASTA_FILE="cython.fasta"
DEV_FASTA_FILE="dev.fasta"

# cd boost.python && python dawg_boost.py > $BOOST_FASTA_FILE
# cd ../
# cd cython       && python dawg_cython.py > $CYTHON_FASTA_FILE
# cd ../
#
# diff $DEV_FASTA_FILE boost.python/$BOOST_FASTA_FILE
# diff $DEV_FASTA_FILE cython/$CYTHON_FASTA_FILE

# Time Cython
cd cython
for (( i = 0; i < 10000; i++ ))
do
	python dawg_cython.py > /dev/null
done
cd ../

# # Time BoostPython
# cd boost.python
# for (( i = 0; i < 10000; i++ ))
# do
# 	python dawg_boost.py > /dev/null
# done
# cd ../

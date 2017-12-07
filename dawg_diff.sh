#!/bin/bash

# Utility script to compare current dawg branch with
# develop outputs
cd ~/repo/dawg/build/debug && git checkout $1 &> /dev/null
cmake -G"Ninja" -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug ../..
ninja -j4

~/repo/dawg/build/debug/src/dawg --seed=1111111 \
    ~/repo/dawg/examples/basic-dna.dawg -o \
    ~/Dropbox/dawg-house/outs/basic-dna_D2.fasta

~/repo/dawg/build/debug/src/dawg --seed=1111111 \
    ~/repo/dawg/examples/basic-dna-2.dawg -o \
    ~/Dropbox/dawg-house/outs/basic-dna-2_D2.fasta

~/repo/dawg/build/debug/src/dawg --seed=1111111 \
    ~/repo/dawg/examples/multiple-models.dawg -o \
    ~/Dropbox/dawg-house/outs/multiple-models_D2.fasta

~/repo/dawg/build/debug/src/dawg --seed=1111111 \
    ~/repo/dawg/examples/recombination.dawg -o \
    ~/Dropbox/dawg-house/outs/recombination_D2.fasta

~/repo/dawg/build/debug/src/dawg --seed=1111111 \
    ~/repo/dawg/examples/segments.dawg -o \
    ~/Dropbox/dawg-house/outs/segments_D2.fasta

#~/repo/dawg/build/debug/src/dawg --seed=1111111 \
#    ~/Dropbox/dawg-house/indels.dawg -o \
#    ~/Dropbox/dawg-house/outs/indels_D2.fasta


diff ~/Dropbox/dawg-house/outs/basicdna_DAWG.fasta ~/Dropbox/dawg-house/outs/basic-dna_D2.fasta
diff ~/Dropbox/dawg-house/outs/segments_DAWG.fasta ~/Dropbox/dawg-house/outs/segments_D2.fasta
diff ~/Dropbox/dawg-house/outs/recomb_DAWG.fasta ~/Dropbox/dawg-house/outs/recombination_D2.fasta
diff ~/Dropbox/dawg-house/outs/multiple-models_DAWG.fasta ~/Dropbox/dawg-house/outs/multiple-models_D2.fasta
#diff ~/Dropbox/dawg-house/outs/indels_DAWG.fasta ~/Dropbox/dawg-house/outs/indels_D2.fasta

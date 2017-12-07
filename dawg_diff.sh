#!/usr/bin/bash

# cd ~/code/dawg

# Simple script that runs dawg examples and validates them
#cmake -G"Ninja" -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug ../..
#ninja -j4

SEED=1111111

build/debug/src/dawg --seed=1111111 \
    examples/basic-dna.dawg -o \
    test_dist/basic-dna_patch.fasta

# ~/code/dawg/build/debug/src/dawg --seed=1111111 \
    # ~/code/dawg/examples/basic-dna-2.dawg -o \
    # ~/Dropbox/dawg-house/outs/basic-dna-2_D2.fasta

build/debug/src/dawg --seed=1111111 \
    examples/multiple-models.dawg -o \
    test_dist/multiple-models_patch.fasta

build/debug/src/dawg --seed=1111111 \
    examples/recombination.dawg -o \
    test_dist/recombination_patch.fasta

build/debug/src/dawg --seed=1111111 \
    examples/segments.dawg -o \
    test_dist/segments_patch.fasta

#~/code/dawg/build/debug/src/dawg --seed=1111111 \
#    ~/Dropbox/dawg-house/indels.dawg -o \
#    ~/Dropbox/dawg-house/outs/indels_D2.fasta


diff ~/code/dawg/test_dist/basicdna.fasta          ~/code/dawg/test_dist/basic-dna_patch.fasta
diff ~/code/dawg/test_dist/segments.fasta          ~/code/dawg/test_dist/segments_patch.fasta
diff ~/code/dawg/test_dist/recomb.fasta            ~/code/dawg/test_dist/recombination_patch.fasta
diff ~/code/dawg/test_dist/multiple-models.fasta   ~/code/dawg/test_dist/multiple-models_patch.fasta
#diff ~/Dropbox/dawg-house/outs/indels_DAWG.fasta ~/Dropbox/dawg-house/outs/indels_D2.fasta

rm ~/code/dawg/test_dist/basic-dna_patch.fasta
rm ~/code/dawg/test_dist/segments_patch.fasta
rm ~/code/dawg/test_dist/recombination_patch.fasta
rm ~/code/dawg/test_dist/multiple-models_patch.fasta

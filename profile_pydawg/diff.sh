#!/bin/bash

# Simple script that diffs the master examples with this branch

SEED=1111111

release/src/dawg --seed=$SEED \
    examples/basic-dna.dawg -o \
    test_dist/basicdna_test.fasta

release/src/dawg --seed=$SEED \
    examples/multiple-models.dawg -o \
    test_dist/multiple-models_test.fasta

release/src/dawg --seed=$SEED \
    examples/recombination.dawg -o \
    test_dist/recombination_test.fasta

release/src/dawg --seed=$SEED \
    examples/segments.dawg -o \
    test_dist/segments_test.fasta

diff -q test_dist/basicdna_develop.fasta                 test_dist/basicdna_test.fasta
diff -q test_dist/segments_develop.fasta                 test_dist/segments_test.fasta
diff -q test_dist/recombination_develop.fasta            test_dist/recombination_test.fasta
diff -q test_dist/multiple-models_develop.fasta          test_dist/multiple-models_test.fasta
#diff ~/Dropbox/dawg-house/outs/indels_DAWG.fasta ~/Dropbox/dawg-house/outs/indels_D2.fasta

rm test_dist/*_test.fasta

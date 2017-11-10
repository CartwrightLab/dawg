#!/bin/bash

# Simple script that diffs the develop examples with this branch
# The develop scripts used this seed
SEED=1111111

release/src/dawg --seed=$SEED \
    examples/basic-dna.dawg -o \
    tests/diff/basicdna_test.fasta

release/src/dawg --seed=$SEED \
    examples/multiple-models.dawg -o \
    tests/diff/multiple-models_test.fasta

release/src/dawg --seed=$SEED \
    examples/recombination.dawg -o \
    tests/diff/recombination_test.fasta

release/src/dawg --seed=$SEED \
    examples/segments.dawg -o \
    tests/diff/segments_test.fasta

diff -q tests/diff/basicdna_develop.fasta        tests/diff/basicdna_test.fasta
diff -q tests/diff/segments_develop.fasta        tests/diff/segments_test.fasta
diff -q tests/diff/recombination_develop.fasta   tests/diff/recombination_test.fasta
diff -q tests/diff/multiple-models_develop.fasta tests/diff/multiple-models_test.fasta

rm tests/diff/*_test.fasta

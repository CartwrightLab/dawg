#!/usr/bin/python3

from os import sys
import PyDawg as pd

print("Hello DAWG, I'm a Python")

# class Bone:
#     __init__(self, bone):
#         self.bone = bone
#
# class Walk:
#     __init__(self, walk):
#         self.walk = walk
#
# class Segments:
#     __init__(self, name, inheritsFrom, *args):
#         self.header   = header
#         self.segments = args
#
# class Trick:
#     __init__(self, name, segments):
#         self.name     = name
#         self.segments = segments
#
# class DawgBuilder:
#     __init__(self, trick, walk, bone):
#         self.trick = trick
#         self.walk  = walk
#         self.bone  = bone

# PyDawg constructor can take in Tricks, Walks, and Bones
akita = pd.pd(
DawgBuilder(
    Trick(name="BasicDnaExample",
        Segments(name="__default__", inheritsFrom="LUCA",
            (tree_tree="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
            root_length=1000, subst_model='hky',
            subst_freqs=[0.2, 0.3, 0.3, 0.2],
            subst_params=[0.2, 1.0], sim_reps=10))),
    Walk(), Bone("fasta:-")))

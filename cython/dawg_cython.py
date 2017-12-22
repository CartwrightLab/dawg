#!/usr/bin/python3
from os import sys
import PyDawg as pd

# dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:/dev/null", 10, 212121)
# dog1.run()
# print("Stats for basic-dna.dawg trick file")
# dog1.trickStats()
# dog1.rand(0, 100))


class TrickSection:

    populateSection(self, name, inheritsFrom, map db):
        self.name = name
        self.inheritsFrom = inheritsFrom
        self.db = db

    compileNewickTree():
        pass

class Trick:

    def buildTrick(self):
        # setup the basic-dna example from dawg without parsing
        ts1 = TrickSection()
        ts1.populateSection(name="basicdna", inheritsFrom="LUCA", [
            "Tree.Tree"="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
            "Subst.Model"="HKY",
            "Subst.Params"="2.0, 1.0",
            "Subst.Freqs"="0.3, 0.2, 0.2, 0.3",
            "Root.Length"="1000",
            "Sim.Reps"="10"
        ])

        ts1.compileNewickTree()
        ts1.


        self.trickSequence = ts1.

if __name__ == '__main__':
    akita = pd.pd()
    akita.trick(name="basicdna",
        [ name="default", inheritsFrom="LUCA",
        "Tree.Tree"="((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);",
        "Subst.Model"="HKY",
        "Subst.Params"="2.0, 1.0",
        "Subst.Freqs"="0.3, 0.2, 0.2, 0.3",
        "Root.Length"="1000",
        "Sim.Reps"="10"
        ])



    print("Hello DAWG, I'm a Python")

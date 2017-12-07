#!/usr/bin/python3
from os import sys
import PyDawg
# import pandas as pd

class TrickSection:

    

# Section commands are redundant
# Global commands are the same for all trick
trickCommands = [
    # Section 1
    TrickSection(subst.model="",
        subst.params="",
        subst.freqs="",
        subst.rate="",),
    TrickSection(),
    TrickSection(),
]


if __name__ == '__main__':
    dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:/dev/null", 10, 212121)
    dog1.run()
    print("Stats for basic-dna.dawg trick file")
    dog1.trickStats()

    numbers = []
    for i in range(1000):
        numbers.append(dog1.rand(0, 100))
    # s = pd.Series(numbers)
    # s.value_counts()
    # df = pd.DataFrame(data=numbers)
    # print(df)
    # pd.DataFrame.hist(data=numbers)

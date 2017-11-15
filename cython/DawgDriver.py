from os import sys
import PyDawg as pd

if __name__ == '__main__':
    dd = pd.PyDawg(b"../examples/basic-dna.dawg", b"fasta:-", 10, 212121)
    dd.run()

#!/usr/bin/python3
from os import sys
import PyDawg
# import pandas as pd

# def createGraph():
#     filename = sys.argv[1]
#     input = pd.read_csv(filename,
#     names=["CompressorName", "CompressorSpeed", \
#         "DecompressSpeed", "Originalsize", "CompressSize", \
#         "ratio", "Filename"], header=None, sep="\s*,\s*")
#
#     averages_speed=input.groupby("CompressorName", as_index=False)["CompressorSpeed"].mean()
#     averages_despeed=input.groupby("CompressorName", as_index=False)["DecompressSpeed"].mean()
#     averages_orgsize=input.groupby("CompressorName", as_index=False)["Originalsize"].mean()
#     averages_comsize=input.groupby("CompressorName", as_index=False)["CompressSize"].mean()
#     averages_ratio=input.groupby("CompressorName", as_index=False)["ratio"].mean()
#
#     table1=averages_speed.drop("CompressorSpeed", axis=1)
#     table2=averages_speed.drop("CompressorName", axis=1)
#     table3=averages_despeed.drop("CompressorName", axis=1)
#     table4=averages_orgsize.drop("CompressorName", axis=1)
#     table5=averages_comsize.drop("CompressorName", axis=1)
#     table6=averages_ratio.drop("CompressorName", axis=1)
#
#     filenamecol = input["Filename"]
#     final = pd.concat([table1, table2, table3, table4, table5, table6], axis=1)
#
#     final_col = pd.concat([final, filenamecol], ignore_index=True, axis=1)
#
#     x = final_col.dropna()
#     y = pd.DataFrame.to_csv(x, header=False, index=False)
#     with pd.option_context('display.max_rows', 60000, 'display.max_columns', 1000, 'display.width', 1000):
#         print(y)

if __name__ == '__main__':
    dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:-", 10, 212121)
    dog2 = PyDawg.PyDawg(212121)
    dog2.rand(0, 100)
    # dog2.printSections()
    # dog2.addSection("SectionA")
    # dog2.printSection("SectionA")
    # section name (the key), inheritance, and a list of values
    # dog2.setSection(b"SectionA", b"LUCA", "Root.Code = 2121")
    alignmentBuffer = PyDawg.trick("").walk("").fetch("") # can pass a treat in fetch to make it go faster

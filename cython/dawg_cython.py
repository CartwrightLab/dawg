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

recombinationStr = """
[[-]]
Subst.Model = jc
Root.Segment = 1
Root.Length = 60
Tree.Tree = ((A:0.02,B:0.02):0.2,(C:0.02):0.2);
[[-]]
Root.Code = 1
Root.Segment = 0
Root.Length = 60
Tree.Tree = ((A:0.02):0.2,(B:0.02,C:0.02):0.2);
[[-]]
Root.Segment = 2
"""

basicDnaStr = """
[Tree]
Tree = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);"
[Subst]
Model = "HKY"
Params = 2.0, 1.0
Freqs = 0.3, 0.2, 0.2, 0.3
[Root]
Length = 1000
[Sim]
Reps = 10
"""

if __name__ == '__main__':
    dog1 = PyDawg.PyDawg(b"../examples/basic-dna.dawg", b"fasta:/dev/null", 10, 212121)
    dog1.run()
    print("Stats for basic-dna.dawg trick file")
    dog1.trickStats()
    # dog1.rand(0, 100)

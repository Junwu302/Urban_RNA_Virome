# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:29:57 2024

@author: Alexander G. Lucaci, Ph.D.
"""

import os
import sys
from glob import glob
import pandas as pd

outputCSV = "merged-phyla-fadeSiteTableWithEBF.csv"
outputCSVfiltered = "filtered-merged-phyla-fadeSiteTableWithEBF.csv"

def mergeCSVs(files: list, outputCSV):
    print(files)
    holder = []
    for _ in files:
        df = pd.read_csv(_)
        df["phylum"] = _.replace(".afa.fasta_fadeSiteTableWithEBF.csv","").replace(".fasta_fadeSiteTableWithEBF.csv", "")
        holder.append(df)
    # end for
    outputDF = pd.concat(holder)
    outputDF.to_csv(outputCSV)
    return outputDF
# end method
        
if not os.path.exists(outputCSV):
    df = mergeCSVs(glob("*.csv"), outputCSV)
else:
    df = pd.read_csv(outputCSV)
    
subsThreshold = 4

filterDF = df[df["subs"] > subsThreshold]

filterDF.reset_index()

filterDF.to_csv(outputCSVfiltered)

# End of file

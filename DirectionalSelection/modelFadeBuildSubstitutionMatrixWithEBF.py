# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:56:20 2024

@author: Alexander G. Lucaci, Ph.D.
"""

# =============================================================================
# Imports
# =============================================================================

import os
import sys
import pandas as pd
from glob import glob
import json
import numpy as np
import re

# =============================================================================
# Declares
# =============================================================================

fadeJsons = glob(os.path.join("..", "*.json"))

outputSuffixCSV = "_fadeSiteTableWithEBF.csv"

# =============================================================================
# Helper functions
# =============================================================================

def processJson(JsonFile: str):
    """
    """
    print("# Processing:", JsonFile)
    _ = open(JsonFile, "r")
    data = json.load(_)
    return data
#end method

"""
1 ['A->I(1)S(3)T(1)', 'E->A(2)D(11)G(5)I(2)K(2)L(1)N(2)P(2)Q(4)R(3)S(3)T(2)V(1)W(1)Y(1)', 'F->A(1)W(1)', 'G->A(3)E(1)N(1)S(3)', 'H->D(2)E(1)F(1)G(4)I(1)K(1)L(2)N(3)Q(2)R(2)S(1)T(1)W(1)', 'I->D(1)F(1)H(1)L(3)M(1)V(7)', 'K->A(1)E(2)M(1)V(2)Y(1)', 'N->D(1)E(3)G(1)K(2)Q(1)R(2)S(1)', 'R->H(1)', 'T->A(1)C(1)D(4)E(3)F(1)G(3)H(1)I(5)K(1)L(5)M(2)N(3)R(5)S(9)V(3)W(2)Y(2)', 'V->A(1)C(1)I(3)Q(1)R(1)T(2)', 'Y->F(2)']
"""

def targetParser(targetData: str):
    # Input: A(3)G(3)N(2)T(4)
    # Output: a list of those [["A", 3], ["G", 3], ["N", "2"], etc]
    matches = re.findall(r'([A-Z])\((\d+)\)', targetData)
    parsed_list = [[letter, int(number)] for letter, number in matches]
    return parsed_list
#end method

def processSites(JsonData: str):
    siteData = JsonData.get("site annotations", np.nan)
    assert siteData != np.nan
    #assert len(siteData["site annotations"]["0"]) = jsonData["input"]["number of sites"]
    print("# Examining:", len(siteData["site annotations"]["0"]), "sites")
    headers = siteData.get("headers", np.nan)
    siteAnnotations = siteData["site annotations"]["0"]
    siteDict = {}
    count = 1
    print("# Looping over sites ...")
    for n, item in enumerate(siteAnnotations):
        #print(n, item)
        composition = item[0]
        annotations = item[1]
        site_detailed_annotations = annotations.split(", ")
        #(n+1, site_detailed_annotations)
        #print("# Looping over detailed site annotations ...")
        for _ in site_detailed_annotations:
            x = _.split("->")
            if len(x) == 1: continue
            #print("x", x)
            source = x[0]
            exchanges = x[1]
            parsedExchanges = targetParser(exchanges)
            #print(parsedExchanges)
            ebf = JsonData["MLE"]["content"][source]["0"][n][3]
            for aa_exchange in parsedExchanges:
                siteDict[count] = {"site": str(n+1), "source": source, "target": aa_exchange[0], "subs": aa_exchange[1], "EBF": ebf}
                count += 1
            # end for
        # end for
    # end for
    print("# Done... returning sites results")
    return siteDict
# end for

# =============================================================================
# MAIN
# =============================================================================
    
for n, file in enumerate(fadeJsons):
    print(f"# Opening json file {file}")
    jsonData = processJson(file)
    print(f"# Processing json data {file}")
    siteDataDict = processSites(jsonData)
    print(f"# Loading dataframe...")
    df = pd.DataFrame.from_dict(siteDataDict, orient='index')
    basename = os.path.basename(file)
    _ = basename.replace(".FADE.json", "") + outputSuffixCSV
    print("# Saving csv file...", _)
    df.to_csv(_, index=False)
    #break
#end for

# =============================================================================
# END OF FILE
# =============================================================================

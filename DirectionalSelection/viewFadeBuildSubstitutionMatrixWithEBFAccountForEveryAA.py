# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:36:43 2024

@author: Alexander G. Lucaci, Ph.D.
"""

# =============================================================================
# Imports
# =============================================================================

import holoviews as hv
import pandas as pd
from holoviews import opts
from random import shuffle
import random
import numpy as np
from glob import glob
import sys
import os

# =============================================================================
# Declares
# =============================================================================

hv.extension('bokeh')

_cmap = "glasbey_dark"

inputCSVs = glob("*_fadeSiteTableWithEBF.csv")

print(f"# Input csv files: {inputCSVs}")

ebf_threshold = 0.0

AminoAcids = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

AA_Holder = []

for aa in AminoAcids:
    _ = pd.DataFrame({"source": aa, "target": aa, "subs": 0}, index = [0])
    AA_Holder.append(_)
# end for

df_AA = pd.concat(AA_Holder)

outputDirectory = "images"

try:
    os.mkdir(outputDirectory)
except:
    pass
# end try

# Custom discrete colormap for the amino acid alphabet
amino_acid_colors = {
    'A': '#1f77b4', 'R': '#ff7f0e', 'N': '#2ca02c', 'D': '#d62728', 'C': '#9467bd',
    'E': '#8c564b', 'Q': '#e377c2', 'G': '#7f7f7f', 'H': '#bcbd22', 'I': '#17becf',
    'L': '#aec7e8', 'K': '#ffbb78', 'M': '#98df8a', 'F': '#ff9896', 'P': '#c5b0d5',
    'S': '#c49c94', 'T': '#f7b6d2', 'W': '#c7c7c7', 'Y': '#dbdb8d', 'V': '#9edae5'
}

colors = [amino_acid_colors[x] for x in amino_acid_colors.keys()]

# =============================================================================
# Helper functions
# =============================================================================

def rotate_label(plot, element):
    white_space = "      "
    angles = plot.handles['text_1_source'].data['angle']
    characters = np.array(plot.handles['text_1_source'].data['text'])
    plot.handles['text_1_source'].data['text'] = np.array([x + white_space if x in characters[np.where((angles < -1.5707963267949) | (angles > 1.5707963267949))] else x for x in plot.handles['text_1_source'].data['text']])
    plot.handles['text_1_source'].data['text'] = np.array([white_space + x if x in characters[np.where((angles > -1.5707963267949) | (angles < 1.5707963267949))] else x for x in plot.handles['text_1_source'].data['text']])
    angles[np.where((angles < -1.5707963267949) | (angles > 1.5707963267949))] += 3.1415926535898
    plot.handles['text_1_glyph'].text_align = "center"
# end method

# =============================================================================
# Load data
# =============================================================================

# =============================================================================
# Preprocessing
# =============================================================================

def preprocess(localdf):
    # Create a Chord object from the DataFrame
    chord = hv.Chord(localdf)
    # Convert DataFrame to Holoviews Dataset
    edges = hv.Dataset(localdf, ['source', 'target'], ['subs'])
    return chord, edges
# end method

# =============================================================================
# Save chart
# =============================================================================

def saveDiagram(chord, file, site, output):
    global outputDirectory
    hv.renderer('bokeh').theme = 'light_minimal'
    print("# Output will be", output)
    hv.save(chord, output, fmt='png', dpi=300)
# end method

# =============================================================================
# Create the chord diagram
# =============================================================================

def createChordDiagram(lcl_df, file, site, output):
    global _cmap
    chord, edges = preprocess(lcl_df)
    #print("# Creating chord")
    # Create the chord diagram
    chord = hv.Chord(edges).opts(
        opts.Chord(
            cmap=_cmap,
            edge_color='source',
            labels='index', 
            #labels='source', 
            node_color='index',
            #node_color='source',
            node_size=20,
            width=900,
            height=900,
            label_text_font_size='18pt',  # Adjust label font size for better visibility
            label_text_color='index',  # Adjust label text color for better contrast
            edge_line_width=hv.dim('subs').norm() * 10,  # Adjust thickness based on 'value'
            show_legend=False,
            edge_alpha=0.7,
            #edge_line_color='black', node_line_width=2
            hooks=[rotate_label]
        )
    )
    chord.options(cmap=['#0000ff', '#8888ff', '#ffffff', '#ff8888', '#ff0000'], colorbar=False, width=400)
    # Save chord
    #print("# Saving file...")
    saveDiagram(chord, file, site, output)
# end method

# =============================================================================
# Main subroutine
# =============================================================================

for file in inputCSVs:
    print(f"# Processing... {file}")
    df = pd.read_csv(file)
    df_statsig = df.copy()
    setOfSites = set(df_statsig["site"])
    setOfSitesList = list(setOfSites)
    for site in setOfSitesList:
        output = file.replace(".csv", "") + "_" + str(site) + ".png"
        output = os.path.join(outputDirectory, output)
        print(f"# Checking for existing file {output}")
        if os.path.exists(output): continue
        print(f"# Creating chart for site {site}")
        tempdf = df[df["site"] == site].copy()
        tempdf.sort_values(by=['site'], ascending=True, inplace=True)
        tempdf.drop(columns=["EBF", "site"], inplace=True)
        tempdf["subs"] = tempdf["subs"] * 100 # ADJUSTMENT!
        #break
        createChordDiagram(tempdf, os.path.basename(file), site, output)
    # end for
    #break
# end for

# =============================================================================
# End of file
# =============================================================================



#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Imports
# =============================================================================

from Bio import SeqIO
import sys
import argparse
from tqdm import tqdm

# =============================================================================
# Declares
# =============================================================================

PROTEIN_FASTA = sys.argv[1]
OUTPUT = sys.argv[2]

# =============================================================================
# Main subroutine.
# =============================================================================

def clean(desc: str):
    name = desc.replace("|", "_")
    name = name.replace("-", "_")
    name = name.replace("=", "_")
    name = name.replace(".", "_")
    name = name.replace(" ", "_")
    name = name.replace("__", "_")
    return name
#end method

results = []

with open(PROTEIN_FASTA, "r") as prot_handle:
    for n, record in tqdm(enumerate(SeqIO.parse(prot_handle, "fasta"))):
        record.description = ""
        record.id = clean(record.id)
        results.append(record)
    #end for
#end with

# Write
SeqIO.write(results, OUTPUT, "fasta")

# =============================================================================
# End of file
# =============================================================================




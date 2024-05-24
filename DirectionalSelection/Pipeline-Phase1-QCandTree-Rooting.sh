#!/bin/bash

files=(data/*.afa)

for key in "${files[@]}"
do
   python3 rename_species.py $key $key.fasta
   bash run_fasttree.sh $key.fasta
done

# Rooting, manually
# This is done on our outgroup for each clade
# We used the website phylotree.hyphy.org

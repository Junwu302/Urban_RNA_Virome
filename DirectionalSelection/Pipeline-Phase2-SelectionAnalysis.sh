#!/bin/bash

# Declares --------------------------------------------------------------------
WD=""

# Software
# NOTE: You may need to modify these variables to be full-paths to the software.
HYPHY="hyphy/bin/hyphy"
HYPHYMPI="hyphy/bin/HYPHYMPI"
RESOURCES="hyphy/res"
NP=8

# HPC Environment-specific ----------------------------------------------------
module load openmpi-4.0.2-gcc-8.2.0-fjrkibe

# Helper functions ------------------------------------------------------------

function run_fade_internal () {
    inputFasta=$1
    inputTree=$2
    outputJson=$inputFasta.FADE.Internal.json
    if [ ! -f $outputJson ]; 
    then
        cmd="mpirun -np $NP $HYPHYMPI LIBPATH=$RESOURCES FADE --alignment $1 --tree $2 --branches Internal --output $outputJson --cache $1.FADE.Internal.cache"
        echo $cmd
        eval $cmd
    fi
}

function run_fade () {
    inputFasta=$1
    inputTree=$2
    outputJson=$inputFasta.FADE.json
    if [ ! -f $outputJson ];
    then
        cmd="mpirun -np $NP $HYPHYMPI LIBPATH=$RESOURCES FADE --alignment $1 --tree $2 --output $outputJson --cache $1.FADE.cache"
        echo $cmd
        eval $cmd
    fi
}

# MAIN ------------------------------------------------------------------------

files=($WD/*.fasta)

for key in "${files[@]}"
do
  echo "FASTA in array is: $key"

  if test -f $key; then
      echo "FASTA exists."
  else
      echo "ERROR! Fasta does not exist"
      #break
  fi
  
  tree="$key.FastTree.rooted.nwk"
  echo "TREE is: $key"

  if test -f $tree; then
      echo "TREE exists."
  else
      echo "ERROR! Tree does not exist"
  fi
 
  echo "Running FADE on all branches..."
  run_fade $key $tree

  echo ""
done


# End of file

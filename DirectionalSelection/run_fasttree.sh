# Description
# Declares
# Software
FASTTREE="FastTree"

run_iqtree () {
    inputFasta=$1
    cmd="$IQTREE -s $inputFasta -nt AUTO"
    echo $cmd
    echo ""
    eval $cmd
}

run_fasttree () {
    inputFasta=$1
    #cmd="$FASTTREE -gamma -nt -gtr $inputFasta > $inputFasta.FastTree.nwk"
    cmd="$FASTTREE -gamma -wag -log logfile $inputFasta > $inputFasta.FastTree.nwk"
    echo $cmd
    echo ""
    eval $cmd
}

run_fasttree $1

# End of file

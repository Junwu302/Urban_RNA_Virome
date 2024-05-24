
# Directional selection analysis of Metasub results

This README describes the steps and software used to conduct directional selection analysis for the MetaSub manuscript.

# Directional selection workflow

## Preliminary steps

Create your working directory and place data in the `data` folder, for this analysis we used all of the available clade data such as 'dupl.afa', 'kitr.afa', 'pisu.afa', etc. The file ending '.afa' is important as we will use this downstream in the analysis.

## Software installation

We rely on the HyPhy (hyphy.org) selection analysis software suite for this analysis. Clone and install the software from the GitHub repository (https://github.com/veg/hyphy) otherwise `conda` installation is also available. We rely on the 'holoviews' python package for visualization you will need to install it into your environment.

## Workflow

The workflow is split into two parts: Phase 1 and Phase 2. 

### Phase 1 - QC and Rooting the Phylogenetic tree

For our analysis, we modify our sequence id's to comply with downstream analysis standard by modifying characters like `|,-,/,=` to underscores. Next, we infer a amino-acid sequence phylogeny using FastTree. The inferred phylogenetic tree is then outrooted manually using the website phylotree.hyphy.org. Run Phase 1 via `bash Pipeline-Phase1-QCandTree-Rooting.sh`

### Phase 2 - Run selection analysis

The rooted phylogenetic tree and amino-acid multiple sequence alignment are used in the FADE (http://vision.hyphy.org/fade) analysis for directional selection. The FADE methods software implementation is available here https://github.com/veg/hyphy/blob/master/res/TemplateBatchFiles/SelectionAnalyses/FADE.bf 

Run Phase 2 via `bash Pipeline-Phase2-SelectionAnalysis.sh`

## Post-processing

Downstream processing of FADE results into visualizations and tables is available. 
 
The python script `modelFadeBuildSubstitutionMatrixWithEBF.py` generates tables from the FADE JSON results for each clade. Call this with `python modelFadeBuildSubstitutionMatrixWithEBF.py`. You will find tables such as: `dupl.afa.fasta_fadeSiteTableWithEBF.csv` and `kitr.fasta_fadeSiteTableWithEBF.csv`

The python script `viewFadeBuildSubstitutionMatrixWithEBFAccountForEveryAA.py` generates chord diagrams for each amino-acid site for each clade. This script will create an 'images' directory which will be populated with the resulting figures. Note, this script takes a while to run. You can call it the same way as the script above.

The python script `summaryStatisticsOnEBFTable.py` is used to merge all of the tables and was used to create our supplementary table.

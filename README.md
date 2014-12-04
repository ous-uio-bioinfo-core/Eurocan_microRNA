Eurocan_microRNA
================

Code used to analyse microRNA expression array data in relation to the findings in Volinia et al.


### parseAFE.r
Script used to parse the Agilent Feature Extraction Files. Creates a object of type uRNAlist to be used later. Only AFE-filse from AHUS are parsed. The UCAM dataset is provided in a processed state.

### read_input.r
Loads the two datasets and sampleannotation, and does some filtering and harmonization. 

### quality_control.rmd
Some quality plotting of the data. Signal distributions, clustring based on data set origin or different grouping labels.

### differential_analysis.rmd
Calculates mean group differences and p-values for groups based on tissue type, pam50 classification or IHC.

### volinia_overlap.rmd
Compares our results with the one replorted in Volinia et al.

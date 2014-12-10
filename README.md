Eurocan_microRNA
================

Code used to analyse microRNA expression array data in relation to the findings in Volinia et al.
The code is split into several smaller files in order to make them more manageable.



### parseAFE.r
Script used to parse the Agilent Feature Extraction Files. Creates a object of type uRNAlist to be used later. Only AFE-files from AHUS are parsed. Will take about 5 minutes to run on a normal computer. The UCAM data set is provided in an already processed state.

### read_input.r
Loads the two data sets and sampleannotation, and does some filtering and harmonization. 

### quality_control.rmd
Some quality control and overview plotting of the data. Signal distributions, clustering based on data set origin or different grouping labels.

### differential_analysis.rmd
Calculates mean group differences and p-values for groups based on tissue type, pam50 classification or IHC. Gene lists for the different comparisons are produced.

### end_result.rmd
Compare our results from the two analysis approaches, meta and merged. Compare our combined results with the one reported in Volinia et al. And finally compare all the reported microRNA against a curated database of microRNA.
 

## NB, not all nececary input files are here!.
Neither the Agilent Feature Extraction Files for the AHUS experiment or the processed matrix from UCAM are in this github-repository. In order to reproduce the reports completely these files need to be obtained and put in a folder named "not_in_github". The same folder is also used to store some intermediate R-objects.
 
 
### Commands to reproduce the reports and results

source("parseAFE.r") # takes about 5 minutes

source("read_input.r") #rest is quick
knit2html("quality_control.rmd")
knit2html("differential_analysis.rmd")
knit2html("end_results.rmd")

parseAFE.r creates a binary file and each script there after can be run individually.


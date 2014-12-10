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
Calculates mean group differences and p-values for groups based on tissue type, pam50 classification or IHC. Genelists for the different comparisons are produced.

### end_result.rmd
Compare our results from the two analysis approaches, meta and merged. Compare our combined results with the one reported in Volinia et al. And finnaly compare alle the reported microRNA against a curated database of microRNA.
 

## NB, not all nececary input files are here!.
Neither the Agilent Feature Extraction Files for the AHUS experiment or the processed matrix from UCAM are in this github-repository. In order to reproduce the reports completly these files need to be obtained and put in a folder named "not_in_github". The same folder is also used to store some intermediate R-objects.
 
 
### Commands to reproduce the reports and results

# source("parseAFE.r") # takes two hours!!
# command to make genelistfiles and html-reports.
# takes little time
source("read_input.r")
knit2html("quality_control.rmd")
knit2html("differential_analysis.rmd")
knit2html("end_results.rmd")



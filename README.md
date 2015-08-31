Eurocan_microRNA
================

This repository consist of R-code used to produce the results reported in TITLE OF OUR PULICATION. 


Two sources of breast cancer microarray microRNA data are analysed, from "AHUS" and "UCAM" (Dvinge et al.), and microRNA that differ between states are reported and compared to the findings in Volinia et al. At last our resulting microRNA are checked towards a curated list of microRNAs reported by Fromm et al.

The code is split into several smaller files in order to make them more manageable.



### parseAFE.r
Script used to parse the Agilent Feature Extraction Files from AHUS. Creates a object of type uRNAlist to be used later. Only AFE-files from AHUS are parsed. Will take about 5 minutes to run on a normal computer. The UCAM data set is provided in an already processed state.

### read_input.r
Loads the two data sets and sampleannotation, and does some filtering and harmonization. 

### quality_control.rmd ( and html)
Some quality control and overview plotting of the data. Signal distributions, clustering based on data set origin or different grouping labels.

### differential_analysis_merged.rmd ( and html)
Calculates mean group differences and p-values for groups based on tissue type, pam50 classification or IHC. Gene lists for the different comparisons are produced. Limma is used and the two data sets are treated as batches and blocked for in Limma.

### differential_analysis_meta.rmd ( and html)
Calculates mean group differences and p-values for groups based on tissue type, pam50 classification or IHC. Gene lists for the different comparisons are produced. Limma is used individually on the two data sets and then a meta analysis is performed based on the p-values.


### end_result.rmd ( and html)
Compares our results from the two analysis approaches, meta and merged. Then compares our combined results with the one reported in Volinia et al. And finally compares all the reported microRNA against a curated database of microRNA.
 

### NB! not all necessary input files are here!.
Neither the Agilent Feature Extraction Files for the AHUS experiment or the processed matrix from UCAM are in this github-repository. In order to reproduce the reports these files need to be obtained and put in a folder named "not_in_github". The same folder is also used to store some intermediate R-objects and old reports.
 
### not_in_github/
- *Agilent_ncRNA_60k_normalised_miRNA_expression_ORIGINAL.txt*, the processed data file made by Dvinge et al. (UCAM) obtainable from ...
- *AHUS_AFE/*, the raw Agilent feature extraction files from AHUS obtainable from ....
 
 
### Pam50 classification not included
The pam50 classifications we use. are based on mRNA microarray data from the same cohorts. Neither the data or script to do this classification are provided here. The pam50 labels are presented and used as any other clinical label. Be aware that the pam50 classification we made for the UCAM data, differ somewhat from the classification in Dvinge et al.


### Commands to reproduce the reports and results
Several libraries might be needed to be installed. 

source("parse_AFE_files.r") # takes about 5 minutes, creates a binary file and each script there after can be run individually.

source("read_input.r") #rest is quick
library(knitr)
knitr::opts_chunk$set(error = FALSE) # halt on error
knit2html("quality_control.rmd")
knit2html("differential_analysis_merged.rmd")
knit2html("differential_analysis_meta.rmd")
knit2html("end_results.rmd")


## References for all scripts

Breast cancer signatures for invasiveness and prognosis defined by deep sequencing of microRNA.  
Volinia S, Galasso M, Sana ME, Wise TF, Palatini J, Huebner K, Croce CM.  
Proc Natl Acad Sci U S A. 2012 Feb 21;109(8):3024-9. doi: 10.1073/pnas.1200010109. Epub 2012 Feb 6.

The shaping and functional consequences of the microRNA landscape in breast cancer.
Dvinge H, Git A, Gräf S, Salmon-Divon M, Curtis C, Sottoriva A, Zhao Y, Hirst M, Armisen J, Miska EA, Chin SF, Provenzano E, Turashvili G, Green A, Ellis I, Aparicio S, Caldas C.
Nature. 2013 May 16;497(7449):378-82. doi: 10.1038/nature12108. Epub 2013 May 5.

A Uniform System For The Annotation Of Human microRNA Genes And The Evolution Of The Human microRNAome. Bastian Fromm, Tyler Billip, Lorenzo Sempere, Morten Johansen, Kjersti Flatmark, Eivind Hovig & Kevin J. Peterson

  
  R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing,
  Vienna, Austria. URL http://www.R-project.org/
  
  Smyth, GK (2005). Limma: linear models for microarray data. In: 'Bioinformatics and Computational Biology Solutions
  using R and Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, pages
  397-420.
  
  Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and
powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B 57,
289–300.

Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.

Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD.(2012) The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics. 2012 Mar 15;28(6):882-3.

  Yihui Xie (2013). knitr: A general-purpose package for dynamic report generation in R. R package version 1.5.

  Yihui Xie (2013) Dynamic Documents with R and knitr. Chapman and Hall/CRC. ISBN 978-1482203530

  Yihui Xie (2013) knitr: A Comprehensive Tool for Reproducible Research in R. In Victoria Stodden, Friedrich Leisch and
  Roger D. Peng, editors, Implementing Reproducible Computational Research. Chapman and Hall/CRC. ISBN 978-1466561595
  
  RStudio Team (2012). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/.





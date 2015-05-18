Ad hoc counting for Bastian
========================================================
2015-05-18 15:58:27


<br/>
<br/>

## Introduction

This file is not part of our Eurocan microRNA project, but special made for Bastian. 

<br/>
<br/>

## Input data

Setting some dependencies

```r
library(limma)
library(xtable)

set.seed(100)
if(!exists("inputisread"))
	source("read_input.r")
```

<br/>
<br/>

Reading Bastian´s list of microRNAs with MIMAT and status.


```r
	curatedfn = "curated_microRNA_MIMAT.csv"
  orgcuratedtab  = read.table(paste("not_in_github", "/", curatedfn, sep=""), sep="\t", 
  												 header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="")
	orgcuratedtab$miRNA = paste("hsa-", orgcuratedtab$miRNA, sep="")
	
	curatedtab = orgcuratedtab[,c(1,2, 4, 3)]
	names(curatedtab)[3]= "MIMAT"
	names(curatedtab)[4]= "MI"
	
	curatedtab = rbind(curatedtab, 
										 #setNames( orgcuratedtab[,c(1,2,4)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,5,3)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,6,3)], names(curatedtab)) )
	curatedtab = curatedtab[curatedtab$MIMAT != "", ]

	table(orgcuratedtab$status)
```

```
## 
##   purple rejected    white   yellow 
##      302      969      484        4
```

The above table is the count of the different types from Bastian´s file.

<br/>
<br/>

## microRNA in Caldas et al.

The array they used was custom made with a lot of *CRI* probes. In total it has **8150** unique microRNA IDs. The breakdown based on type (as deduced by the 3 first letters:)

```r
table(substring(ucam_genes, 1, 3))
```

```
## 
##  bkv  Bla  CRI  dmr  ebv  hbv  hcm  hiv  hsa  hsv  hur  jcv  ksh  miR  mr_ 
##    2    1 6902    6   44    3   17    4 1083   48    5    2   24    1    1 
##  NC1  NC2  Neg 
##    2    4    1
```
<br/>
<br/>

Mapping names from the microarray to Bastians list via MIMATID.


```r
comb = data.frame(uniqueID=union(unfilteredUCAMMIMAT, unfilteredAHUSMIMAT), stringsAsFactors=FALSE)
comb$lookupID = unlist(lapply(comb$uniqueID, FUN=function(x)strsplit(x, "_")[[1]][1]))
comb$MI = curatedtab[match(comb$lookupID, curatedtab$MIMAT), "MI"]
comb$status = curatedtab[match(comb$lookupID, curatedtab$MIMAT), "status"]
comb$name = curatedtab[match(comb$lookupID, curatedtab$MIMAT), "miRNA"]
comb$AHUS = "notonarray"
comb$AHUS[(comb$uniqueID %in% unfilteredAHUSMIMAT)]="notexpressed"
comb$AHUS[(comb$uniqueID %in% filteredAHUSMIMAT)]="expressed"
comb$UCAM = "notonarray"
comb$UCAM[(comb$uniqueID %in% unfilteredUCAMMIMAT)]="notexpressed"
comb$UCAM[(comb$uniqueID %in% filteredUCAMMIMAT)]="expressed"
comb = comb[!is.na(comb$status),]
#comb[1:10,]
#table(comb[comb$status=="purple", c("AHUS", "UCAM")])


tab = rbind(nonexpressed=table(comb$status[comb$UCAM=="notexpressed"]),
						expressed=table(comb$status[comb$UCAM=="expressed"]))

tab = cbind(tab, sum=rowSums(tab))
perctab = round(tab/rowSums(tab) * 100 * 2, 0)
x = paste(tab,paste("(", perctab, "%)", sep=""), sep=" ")
nicetab = matrix(x, nrow=nrow(perctab), ncol=ncol(perctab), dimnames=dimnames(perctab))

print( xtable(nicetab, caption="", digits=0), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> purple </th> <th> rejected </th> <th> white </th> <th> yellow </th> <th> sum </th>  </tr>
  <tr> <td align="right"> nonexpressed </td> <td> 81 (18%) </td> <td> 174 (39%) </td> <td> 190 (43%) </td> <td> 1 (0%) </td> <td> 446 (100%) </td> </tr>
  <tr> <td align="right"> expressed </td> <td> 40 (7%) </td> <td> 106 (18%) </td> <td> 439 (75%) </td> <td> 2 (0%) </td> <td> 587 (100%) </td> </tr>
   </table>
The above table might count a microRNA two times. Several different microRNA names from the array might end up link to the same microRNA listed in Bastian´s file. Both arms of the miroRNA can be printed, and sometimes one is expressed, while the other is not. This is difficult to count in an easy understandable way. If I counted unique microRNA´s, some will be both expressed and non-expressed.
<br/>
<br/>
The above numbers differed from what we got from the combined data set, especially for the purple type. Therefore, I do the same table for the AHUS only data:


```r
tab = rbind(nonexpressed=table(comb$status[comb$AHUS=="notexpressed"]),
						expressed=table(comb$status[comb$AHUS=="expressed"]))
tab = cbind(tab, sum=rowSums(tab))
perctab = round(tab/rowSums(tab) * 100 * 2, 0)
x = paste(tab,paste("(", perctab, "%)", sep=""), sep=" ")
nicetab = matrix(x, nrow=nrow(perctab), ncol=ncol(perctab), dimnames=dimnames(perctab))

print( xtable(nicetab, caption="", digits=0), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> purple </th> <th> rejected </th> <th> white </th> <th> yellow </th> <th> sum </th>  </tr>
  <tr> <td align="right"> nonexpressed </td> <td> 70 (13%) </td> <td> 112 (21%) </td> <td> 354 (66%) </td> <td> 2 (0%) </td> <td> 538 (100%) </td> </tr>
  <tr> <td align="right"> expressed </td> <td> 3 (1%) </td> <td> 35 (13%) </td> <td> 236 (86%) </td> <td> 1 (0%) </td> <td> 275 (100%) </td> </tr>
   </table>

```r
write.table(comb, file="temp/microRNAwithstatus.txt",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="")

print( paste("Unique MI: ", length(unique(comb$MI))))
```

[1] "Unique MI:  794"

```r
# comb = data.frame(uniqueID=unfilteredcommonMIMAT)
# comb$lookupID = unlist(lapply(unfilteredcommonMIMAT, FUN=function(x)strsplit(x, "_")[[1]][1]))
# comb$status = curatedtab[match(comb$lookupID, curatedtab$MIMAT), "status"]
# comb$name = curatedtab[match(comb$lookupID, curatedtab$MIMAT), "miRNA"]
# comb$AHUS = "unknown"
# comb$AHUS[(comb$uniqueID %in% unfilteredAHUSMIMAT)]="notexpressed"
# comb$AHUS[(comb$uniqueID %in% filteredAHUSMIMAT)]="expressed"
# comb$UCAM = "unknown"
# comb$UCAM[(comb$uniqueID %in% unfilteredUCAMMIMAT)]="notexpressed"
# comb$UCAM[(comb$uniqueID %in% filteredUCAMMIMAT)]="expressed"
# comb=comb[!is.na(comb$status),]
#table(comb[comb$status=="purple", c("AHUS", "UCAM")])
```

The two data sets differ quite a bit in what types of microRNA that are filtered out. It might be an error somewhere in the counting, but I was unable to find it. From the numbers it is also clear that the filtering of the two data sets differ. More microRNAs are lost in the AHUS filtering, and this might explain the loss of almost all the purples.







## References

The shaping and functional consequences of the microRNA landscape in breast cancer.
Dvinge H, Git A, Gräf S, Salmon-Divon M, Curtis C, Sottoriva A, Zhao Y, Hirst M, Armisen J, Miska EA, Chin SF, Provenzano E, Turashvili G, Green A, Ellis I, Aparicio S, Caldas C.
Nature. 2013 May 16;497(7449):378-82. doi: 10.1038/nature12108. Epub 2013 May 5.

  
  R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing,
  Vienna, Austria. URL http://www.R-project.org/

  Yihui Xie (2013). knitr: A general-purpose package for dynamic report generation in R. R package version 1.5.

  Yihui Xie (2013) Dynamic Documents with R and knitr. Chapman and Hall/CRC. ISBN 978-1482203530

  Yihui Xie (2013) knitr: A Comprehensive Tool for Reproducible Research in R. In Victoria Stodden, Friedrich Leisch and
  Roger D. Peng, editors, Implementing Reproducible Computational Research. Chapman and Hall/CRC. ISBN 978-1466561595
  
  RStudio Team (2012). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/.



```r
sessionInfo()
```

```
R version 3.1.1 (2014-07-10)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] MAMA_2.2.1            GeneMeta_1.36.0       gtools_3.4.1         
 [4] multtest_2.20.0       metaMA_2.1            SMVar_1.3.3          
 [7] genefilter_1.46.1     sva_3.10.0            mgcv_1.8-4           
[10] nlme_3.1-120          corpcor_1.6.7         xtable_1.7-4         
[13] RColorBrewer_1.1-2    knitr_1.9             data.table_1.9.2     
[16] plyr_1.8.1            AgiMicroRna_2.14.0    affycoretools_1.36.1 
[19] GO.db_2.14.0          RSQLite_1.0.0         DBI_0.3.1            
[22] AnnotationDbi_1.26.1  GenomeInfoDb_1.0.2    preprocessCore_1.26.1
[25] affy_1.42.3           limma_3.20.9          Biobase_2.24.0       
[28] BiocGenerics_0.10.0  

loaded via a namespace (and not attached):
 [1] acepack_1.3-3.3           affyio_1.32.0            
 [3] annaffy_1.36.0            annotate_1.42.1          
 [5] AnnotationForge_1.6.1     base64enc_0.1-2          
 [7] BatchJobs_1.5             BBmisc_1.9               
 [9] BiocInstaller_1.14.3      BiocParallel_0.6.1       
[11] biomaRt_2.20.0            Biostrings_2.32.1        
[13] biovizBase_1.12.3         bit_1.1-12               
[15] bitops_1.0-6              brew_1.0-6               
[17] BSgenome_1.32.0           Category_2.30.0          
[19] caTools_1.17.1            checkmate_1.5.1          
[21] cluster_2.0.1             codetools_0.2-10         
[23] colorspace_1.2-4          DESeq2_1.4.5             
[25] dichromat_2.0-0           digest_0.6.8             
[27] edgeR_3.6.8               evaluate_0.5.5           
[29] fail_1.2                  ff_2.2-13                
[31] foreach_1.4.2             foreign_0.8-63           
[33] formatR_1.0               Formula_1.2-0            
[35] gcrma_2.36.0              gdata_2.13.3             
[37] geneplotter_1.42.0        GenomicAlignments_1.0.6  
[39] GenomicFeatures_1.16.3    GenomicRanges_1.16.4     
[41] ggbio_1.12.10             ggplot2_1.0.0            
[43] GOstats_2.30.0            gplots_2.16.0            
[45] graph_1.42.0              gridExtra_0.9.1          
[47] GSEABase_1.26.0           gtable_0.1.2             
[49] Hmisc_3.15-0              hwriter_1.3.2            
[51] IRanges_1.22.10           iterators_1.0.7          
[53] KernSmooth_2.23-14        lattice_0.20-30          
[55] latticeExtra_0.6-26       locfit_1.5-9.1           
[57] markdown_0.7.4            MASS_7.3-39              
[59] Matrix_1.1-5              MergeMaid_2.36.0         
[61] metaArray_1.42.0          mime_0.2                 
[63] munsell_0.4.2             nnet_7.3-9               
[65] oligoClasses_1.26.0       PFAM.db_2.14.0           
[67] proto_0.3-10              R.methodsS3_1.7.0        
[69] R.oo_1.19.0               R.utils_1.34.0           
[71] R2HTML_2.3.1              RBGL_1.40.1              
[73] Rcpp_0.11.4               RcppArmadillo_0.4.650.1.1
[75] RCurl_1.95-4.5            ReportingTools_2.4.0     
[77] reshape2_1.4.1            rpart_4.1-9              
[79] Rsamtools_1.16.1          rtracklayer_1.24.2       
[81] scales_0.2.4              sendmailR_1.2-1          
[83] splines_3.1.1             stats4_3.1.1             
[85] stringr_0.6.2             survival_2.38-1          
[87] tools_3.1.1               VariantAnnotation_1.10.5 
[89] XML_3.98-1.1              XVector_0.4.0            
[91] zlibbioc_1.10.0          
```

generation ended 2015-05-18 15:58:28. Time spent 0 minutes .

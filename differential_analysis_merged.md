Finding differentially expressed microRNA in the  AHUS and UCAM data sets
========================================================
2015-06-03 17:44:21


<br/>
<br/>

## Introduction

The aim of the analysis performed in this report is to validate some of the findings from [Volinia et al.](http://www.pnas.org/content/early/2012/02/01/1200010109), 
based on two independent data sets using 3 different Agilent microRNA microarrays from two providers called **AHUS** and  **UCAM**.

More specifically:

- Combine the data into one big sample vs. microRNA matrix
- Do a "find differentially expressed genes"- test between the tissue types DCIS-normal and invasive-DCIS
- With a significance cut-off, report the overlap with what is given in Volinia et al. [(Suppl. table S1 and S2)](http://www.pnas.org/content/suppl/2012/02/02/1200010109.DCSupplemental/pnas.201200010SI.pdf) 
- Do an additional  "find differentially expressed genes"- test between the tissue types DCIS and the individual pam50 subtypes based on classification performed on mRNA samples.

This report is made with RStudio and knitr in order to tie the description, code, plots and results together, and hopefully to limit confusion and adhere to a more reproducible research style. 


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

outputdir = paste( "output-analysis", sep="")
if(!file.exists(outputdir))
	dir.create(outputdir)
difflistdir = paste(outputdir, "/difflists_merged", sep="")
if(!file.exists(difflistdir))
	dir.create(difflistdir)
postfix = ".txt"
```


The datamatrix has 1573 samples and 266 microRNAs.



```r
print(xtable(table(sampleannotation[, c("provider", "tissue_type")]), 
             caption="", digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun  3 17:44:21 2015 -->
<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> benign </th> <th> DCIS </th> <th> invasive </th> <th> normal </th>  </tr>
  <tr> <td align="right"> AHUS </td> <td align="right">   23 </td> <td align="right">    8 </td> <td align="right">   55 </td> <td align="right">   70 </td> </tr>
  <tr> <td align="right"> UCAM </td> <td align="right">    8 </td> <td align="right">   10 </td> <td align="right"> 1283 </td> <td align="right">  116 </td> </tr>
   </table>

The tissue is not evenly distributed among the providers. "benign" is not part of Volinia et al.s results, but it is kept for now since it is useful in assessing batch differences. 

<br/>
<br/>

## Finding differetially expressed microRNAs

Limma is used and both batch and group is included in the analysis.

DCIS vs normal and invasive vs. DCIS comparisons.

```r
# setting output folder

#Limma blocked batch and significance test
group = factor(sampleannotation$tissue_type)
batch = factor(sampleannotation$provider)
 
#  a = batch=="AHUS" & group %in% c("DCIS", "invasive")
#  group[a] = sample(group[a])
#  b = batch=="UCAM" & group %in% c("DCIS", "invasive")
#  group[b] = sample(group[b])

design = model.matrix(~0+group+batch)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix, design)
fit$genes$name = unfilteredcommonname[rownames(common_matrix)]
cont.matrix = makeContrasts ( contrasts=c("DCIS-normal", "invasive-DCIS"), levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DCIS_vs_normal_test = topTable(fit2, coef="DCIS-normal", adjust="BH", number=999)
print( paste("Diff genes with fdr<0.05 for DCIS vs. normal: ",
               sum(DCIS_vs_normal_test$adj.P.Val < 0.05), sep=""))
```

[1] "Diff genes with fdr<0.05 for DCIS vs. normal: 86"

```r
DCIS_vs_normal_difflistfile = paste(difflistdir, 
                                    "/difflist-DCIS-vs-normal",
                                    postfix, sep="")
write.table(DCIS_vs_normal_test, 
            file=DCIS_vs_normal_difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
print(xtable(head(DCIS_vs_normal_test), caption="Top of the result from the DCIS vs normal test", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> Top of the result from the DCIS vs normal test </caption>
<tr> <th> name </th> <th> logFC </th> <th> AveExpr </th> <th> t </th> <th> P.Value </th> <th> adj.P.Val </th> <th> B </th>  </tr>
  <tr> <td> hsa-miR-1305 </td> <td align="right"> 0.681 </td> <td align="right"> 6.564 </td> <td align="right"> 13.003 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 72.627 </td> </tr>
  <tr> <td> hsa-miR-1288 </td> <td align="right"> 0.342 </td> <td align="right"> 6.256 </td> <td align="right"> 12.176 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 63.205 </td> </tr>
  <tr> <td> hsa-miR-21 </td> <td align="right"> 2.407 </td> <td align="right"> 14.925 </td> <td align="right"> 10.902 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 49.709 </td> </tr>
  <tr> <td> hsa-miR-139-5p </td> <td align="right"> -1.276 </td> <td align="right"> 7.033 </td> <td align="right"> -9.693 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 38.083 </td> </tr>
  <tr> <td> hsa-miR-1274b </td> <td align="right"> 1.600 </td> <td align="right"> 11.862 </td> <td align="right"> 9.690 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 38.057 </td> </tr>
  <tr> <td> hsa-miR-486-5p </td> <td align="right"> -1.820 </td> <td align="right"> 7.517 </td> <td align="right"> -7.843 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 22.678 </td> </tr>
   </table>

Many microRNAs seems to differ between DCIS and normal.

Now, the invasive vs. DCIS comparison.

```r
invasive_vs_DCIS_test =  topTable(fit2, coef="invasive-DCIS", adjust="BH", number=999)
  print( paste("Diff genes with fdr<0.05 for invasive vs. DCIS: ",
               sum(invasive_vs_DCIS_test$adj.P.Val < 0.05), sep=""))
```

[1] "Diff genes with fdr<0.05 for invasive vs. DCIS: 10"

```r
invasive_vs_DCIS_difflistfile=paste(difflistdir, 
                                    "/difflist-invasive-vs-DCIS",
                                    postfix, sep="") 
write.table(invasive_vs_DCIS_test, 
            file=invasive_vs_DCIS_difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
print(xtable(head(invasive_vs_DCIS_test), caption="Top of the result from the invasive vs DCIS test", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> Top of the result from the invasive vs DCIS test </caption>
<tr> <th> name </th> <th> logFC </th> <th> AveExpr </th> <th> t </th> <th> P.Value </th> <th> adj.P.Val </th> <th> B </th>  </tr>
  <tr> <td> hsa-miR-1288 </td> <td align="right"> -0.302 </td> <td align="right"> 6.256 </td> <td align="right"> -11.060 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 51.376 </td> </tr>
  <tr> <td> hsa-miR-1305 </td> <td align="right"> -0.546 </td> <td align="right"> 6.564 </td> <td align="right"> -10.704 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 47.793 </td> </tr>
  <tr> <td> hsa-miR-1274b </td> <td align="right"> -0.931 </td> <td align="right"> 11.862 </td> <td align="right"> -5.791 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 9.260 </td> </tr>
  <tr> <td> hsa-miR-15b </td> <td align="right"> 0.715 </td> <td align="right"> 12.638 </td> <td align="right"> 4.750 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 3.893 </td> </tr>
  <tr> <td> hsa-miR-720 </td> <td align="right"> -0.793 </td> <td align="right"> 14.486 </td> <td align="right"> -4.177 </td> <td align="right"> 0.000 </td> <td align="right"> 0.002 </td> <td align="right"> 1.376 </td> </tr>
  <tr> <td> hsa-miR-107 </td> <td align="right"> 0.588 </td> <td align="right"> 10.031 </td> <td align="right"> 3.963 </td> <td align="right"> 0.000 </td> <td align="right"> 0.003 </td> <td align="right"> 0.519 </td> </tr>
   </table>



## Stratification based on pam50 classification.

Pam50 classification has been done on the mRNA samples. Many of the microRNA samples also have a mRNA sample and thus a pam50 classification. 


```r
table(sampleannotation[,c("tissue_type", "pam50_est")], useNA="ifany")
```

```
##            pam50_est
## tissue_type Basallike DCIS Her2 LumA LumB Normallike unknown
##    benign           2    0    0    3    0         12      14
##    DCIS             0   18    0    0    0          0       0
##    invasive       198    0  167  476  385         97      15
##    normal           0    0    0    9    1        105      71
```

We have too few samples for looking at the DCIS to invasive change for each of the pam50 classes. We will treat all DCIS as "DCIS".


```r
pam50 = sampleannotation$pam50_est
pam50[sampleannotation$tissue_type=="DCIS"]="DCIS"
pam50labels = unique(pam50)
pam50labels = pam50labels[!pam50labels %in% c("unknown", "DCIS")]
```


The diff test,

```r
pam50contrasts = paste(pam50labels, "-DCIS", sep="")
to_be_used = sampleannotation$tissue_type %in% c("DCIS", "invasive") & (pam50 != "unknown") # do not use samples without clinical info!

group = factor(pam50[to_be_used])
batch = factor(sampleannotation$provider[to_be_used])

# simple permutation sanity check
#a = batch=="AHUS" & group %in% c("DCIS", pam50labels)
#group[a] = sample(group[a])
#b = batch=="UCAM" & group %in% c("DCIS", pam50labels)
#group[b] = sample(group[b])

design = model.matrix(~0+group+batch)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix[,to_be_used], design)
fit$genes$name = unfilteredcommonname[rownames(common_matrix)]
cont.matrix = makeContrasts ( contrasts=pam50contrasts, levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

for(thiscontrast in pam50contrasts)
{
  thisres = topTable(fit2, coef=thiscontrast, adjust="BH", number=999)

  thisname= gsub("-", "-vs-", thiscontrast)
  thislabel = strsplit(thiscontrast, "-")[[1]][1]

  difflistfile=paste(difflistdir, "/difflist-pam50-",thisname, postfix, sep="")
  write.table(thisres, 
            file=difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  print( paste("Diff genes with fdr<0.05 for ", thislabel, " vs DCIS: ",
               sum(thisres$adj.P.Val < 0.05), sep=""))
  
}
```

```
## [1] "Diff genes with fdr<0.05 for Normallike vs DCIS: 34"
## [1] "Diff genes with fdr<0.05 for LumA vs DCIS: 31"
## [1] "Diff genes with fdr<0.05 for LumB vs DCIS: 43"
## [1] "Diff genes with fdr<0.05 for Her2 vs DCIS: 16"
## [1] "Diff genes with fdr<0.05 for Basallike vs DCIS: 61"
```

This is diff genes between all the DCIS (which is a small and diverse group) and each of the pam50 classes for the invasive samples.

<br/>
<br/>


## Stratification based on Immunohistochemistry (IHC) classification.

Many of the samples are also diagnosed based on immunohistochemistry of HER2, PGR and ER proteins. Four diagnoses are used, 
"HER2neg/ERneg/PGRneg", "HER2neg/ERpos", "HER2pos/ERneg", "HER2pos/ERpos". 
We would like too look for microRNA that are changed between all the DCIS and the different IHC diagnoses of invasive samples.


```r
IHC = sampleannotation$IHC
IHC[sampleannotation$tissue_type=="DCIS"] = "DCIS"
table(sampleannotation[,"tissue_type"], IHC, useNA="ifany")
```

```
##           IHC
##            DCIS HER2neg_ERneg_PGRneg HER2neg_ERpos HER2pos_ERneg
##   benign      0                    5             2             0
##   DCIS       18                    0             0             0
##   invasive    0                  199           940            87
##   normal      0                    0             0             0
##           IHC
##            HER2pos_ERpos unknown
##   benign               0      24
##   DCIS                 0       0
##   invasive            93      19
##   normal               0     186
```


The diff test,

```r
IHClabels = unique(IHC[!IHC %in% c("", "unknown", "DCIS")])
IHCcontrasts = paste(IHClabels, "-DCIS", sep="")
to_be_used =  sampleannotation$tissue_type %in% c("DCIS", "invasive") & (IHC != "unknown") # do not use samples without clinical info!

group = factor(IHC[to_be_used])
batch = factor(sampleannotation$provider[to_be_used])

#simple permutation sanity check
# a = batch=="AHUS" & group %in% c("DCIS", IHClabels)
# group[a] = sample(group[a])
# b = batch=="UCAM" & group %in% c("DCIS", IHClabels)
# group[b] = sample(group[b])

design = model.matrix(~0+group+batch)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix[,to_be_used], design)
fit$genes$name = unfilteredcommonname[rownames(common_matrix)]
cont.matrix = makeContrasts ( contrasts=IHCcontrasts, levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

for(thiscontrast in IHCcontrasts)
{
  thisres = topTable(fit2, coef=thiscontrast, adjust="BH", number=999)

  thisname= gsub("-", "-vs-", thiscontrast)
  thislabel = strsplit(thiscontrast, "-")[[1]][1]

  difflistfile=paste(difflistdir, "/difflist-IHC-",thisname, postfix, sep="")
  write.table(thisres, 
            file=difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  print( paste("Diff genes with fdr<0.05 for ", thislabel, " vs DCIS: ",
               sum(thisres$adj.P.Val < 0.05), sep=""))
  
}
```

```
## [1] "Diff genes with fdr<0.05 for HER2neg_ERpos vs DCIS: 11"
## [1] "Diff genes with fdr<0.05 for HER2pos_ERneg vs DCIS: 19"
## [1] "Diff genes with fdr<0.05 for HER2pos_ERpos vs DCIS: 5"
## [1] "Diff genes with fdr<0.05 for HER2neg_ERneg_PGRneg vs DCIS: 30"
```



**Important note regarding Limma results**
For the diff-tests I ran a few sanity checks where I permuted the group labels whitin batches. For every 3 of 10 run I got a few genes, and once in a while I did get many more. This indicates that some assumtions of the model might be broken to a degree that the resulting p-values should be considered somewhat optimistic.



<br>
<br>


```r
sessionInfo()
```

```
R version 3.1.1 (2014-07-10)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] sva_3.10.0            mgcv_1.8-4            nlme_3.1-120         
 [4] corpcor_1.6.7         xtable_1.7-4          RColorBrewer_1.1-2   
 [7] knitr_1.9             AgiMicroRna_2.14.0    affycoretools_1.36.1 
[10] GO.db_2.14.0          RSQLite_1.0.0         DBI_0.3.1            
[13] AnnotationDbi_1.26.1  GenomeInfoDb_1.0.2    preprocessCore_1.26.1
[16] affy_1.42.3           limma_3.20.9          Biobase_2.24.0       
[19] BiocGenerics_0.10.0   data.table_1.9.2      plyr_1.8.1           

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
[37] genefilter_1.46.1         geneplotter_1.42.0       
[39] GenomicAlignments_1.0.6   GenomicFeatures_1.16.3   
[41] GenomicRanges_1.16.4      ggbio_1.12.10            
[43] ggplot2_1.0.0             GOstats_2.30.0           
[45] gplots_2.16.0             graph_1.42.0             
[47] grid_3.1.1                gridExtra_0.9.1          
[49] GSEABase_1.26.0           gtable_0.1.2             
[51] gtools_3.4.1              Hmisc_3.15-0             
[53] hwriter_1.3.2             IRanges_1.22.10          
[55] iterators_1.0.7           KernSmooth_2.23-14       
[57] lattice_0.20-30           latticeExtra_0.6-26      
[59] locfit_1.5-9.1            markdown_0.7.4           
[61] MASS_7.3-39               Matrix_1.1-5             
[63] mime_0.2                  munsell_0.4.2            
[65] nnet_7.3-9                oligoClasses_1.26.0      
[67] PFAM.db_2.14.0            proto_0.3-10             
[69] R.methodsS3_1.7.0         R.oo_1.19.0              
[71] R.utils_1.34.0            R2HTML_2.3.1             
[73] RBGL_1.40.1               Rcpp_0.11.4              
[75] RcppArmadillo_0.4.650.1.1 RCurl_1.95-4.5           
[77] ReportingTools_2.4.0      reshape2_1.4.1           
[79] rpart_4.1-9               Rsamtools_1.16.1         
[81] rtracklayer_1.24.2        scales_0.2.4             
[83] sendmailR_1.2-1           splines_3.1.1            
[85] stats4_3.1.1              stringr_0.6.2            
[87] survival_2.38-1           tools_3.1.1              
[89] VariantAnnotation_1.10.5  XML_3.98-1.1             
[91] XVector_0.4.0             zlibbioc_1.10.0          
```

generation ended 2015-06-03 17:44:22. Time spent 0 minutes .

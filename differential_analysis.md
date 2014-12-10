Finding differentially expressed microRNA in the  AHUS and UCAM data sets
========================================================



2014-12-10 16:24:46  
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
	{dir.create(outputdir)}else{
		# move old results and add a datastring
		datestr = strsplit( as.character(file.info(outputdir)$ctime), " ")[[1]][1]
    backupname = paste("not_in_github/", outputdir, "_", datestr, sep="")
    x=file.rename(outputdir,backupname)
		dir.create(outputdir)
	}
```

```
## Warning in file.rename(outputdir, backupname): cannot rename file
## 'output-analysis' to 'not_in_github/output-analysis_2014-12-08', reason
## 'Directory not empty'
```

```
## Warning in dir.create(outputdir): 'output-analysis' already exists
```

```r
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
<!-- Wed Dec 10 16:24:46 2014 -->
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
cont.matrix = makeContrasts ( contrasts=c("DCIS-normal", "invasive-DCIS"), levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

DCIS_vs_normal_test = topTable(fit2, coef="DCIS-normal", adjust="BH", number=999)
print( paste("Diff genes with fdr<0.05 for DCIS vs. normal: ",
               sum(DCIS_vs_normal_test$adj.P.Val < 0.05), sep=""))
```

[1] "Diff genes with fdr<0.05 for DCIS vs. normal: 86"

```r
DCIS_vs_normal_difflistfile = paste(outputdir, 
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
<tr> <th> logFC </th> <th> AveExpr </th> <th> t </th> <th> P.Value </th> <th> adj.P.Val </th> <th> B </th>  </tr>
  <tr> <td align="right"> 0.681 </td> <td align="right"> 6.564 </td> <td align="right"> 13.003 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 72.627 </td> </tr>
  <tr> <td align="right"> 0.342 </td> <td align="right"> 6.256 </td> <td align="right"> 12.176 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 63.205 </td> </tr>
  <tr> <td align="right"> 2.407 </td> <td align="right"> 14.925 </td> <td align="right"> 10.902 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 49.709 </td> </tr>
  <tr> <td align="right"> -1.276 </td> <td align="right"> 7.033 </td> <td align="right"> -9.693 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 38.083 </td> </tr>
  <tr> <td align="right"> 1.600 </td> <td align="right"> 11.862 </td> <td align="right"> 9.690 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 38.057 </td> </tr>
  <tr> <td align="right"> -1.820 </td> <td align="right"> 7.517 </td> <td align="right"> -7.843 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 22.678 </td> </tr>
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
invasive_vs_DCIS_difflistfile=paste(outputdir, 
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
<tr> <th> logFC </th> <th> AveExpr </th> <th> t </th> <th> P.Value </th> <th> adj.P.Val </th> <th> B </th>  </tr>
  <tr> <td align="right"> -0.302 </td> <td align="right"> 6.256 </td> <td align="right"> -11.060 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 51.376 </td> </tr>
  <tr> <td align="right"> -0.546 </td> <td align="right"> 6.564 </td> <td align="right"> -10.704 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 47.793 </td> </tr>
  <tr> <td align="right"> -0.931 </td> <td align="right"> 11.862 </td> <td align="right"> -5.791 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 9.260 </td> </tr>
  <tr> <td align="right"> 0.715 </td> <td align="right"> 12.638 </td> <td align="right"> 4.750 </td> <td align="right"> 0.000 </td> <td align="right"> 0.000 </td> <td align="right"> 3.893 </td> </tr>
  <tr> <td align="right"> -0.793 </td> <td align="right"> 14.486 </td> <td align="right"> -4.177 </td> <td align="right"> 0.000 </td> <td align="right"> 0.002 </td> <td align="right"> 1.376 </td> </tr>
  <tr> <td align="right"> 0.588 </td> <td align="right"> 10.031 </td> <td align="right"> 3.963 </td> <td align="right"> 0.000 </td> <td align="right"> 0.003 </td> <td align="right"> 0.519 </td> </tr>
   </table>



## Stratification based on pam50 classification.

Pam50 classification has been done on the mRNA samples. Many of the microRNA samples also have a mRNA sample and thus a pam50 classification. 
We have too few samples for looking at the DCIS to invasive change for each of the pam50 classes. For now treat all DCIS as "DCIS".


```r
pam50labels = unique(sampleannotation$pam50_est)
pam50labels = pam50labels[!pam50labels %in% c("unknown", "DCIS")]
sampleannotation$pam50_est[sampleannotation$tissue_type=="DCIS"]="DCIS"
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


The diff test,

```r
#signchangetab = matrix( nrow=nrow(common_matrix_batchnorm), ncol=length(pam50labels),
#                       dimnames=list(rownames(common_matrix_batchnorm), pam50labels))
dcis_invasive = sampleannotation$tissue_type %in% c("DCIS", "invasive")
signchangetab = as.data.frame(matrix( nrow=nrow(common_matrix), ncol=length(pam50labels)+2), stringsAsFactors=FALSE)
colnames(signchangetab) = c("microRNA", "MIMATID", pam50labels)
rownames(signchangetab) = rownames(common_matrix)
signchangetab[,"MIMATID"] = rownames(common_matrix)

pam50contrasts = paste(pam50labels, "-DCIS", sep="")


group = factor(sampleannotation$pam50_est[dcis_invasive])
batch = factor(sampleannotation$provider[dcis_invasive])

# simple permutation sanity check
#a = batch=="AHUS" & group %in% c("DCIS", pam50labels)
#group[a] = sample(group[a])
#b = batch=="UCAM" & group %in% c("DCIS", pam50labels)
#group[b] = sample(group[b])

design = model.matrix(~0+group+batch)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix[,dcis_invasive], design)
cont.matrix = makeContrasts ( contrasts=pam50contrasts, levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

for(thiscontrast in pam50contrasts)
{
  thisres = topTable(fit2, coef=thiscontrast, adjust="BH", number=999)

  thisname= gsub("-", "-vs-", thiscontrast)
  thislabel = strsplit(thiscontrast, "-")[[1]][1]

  difflistfile=paste(outputdir, "/difflist-pam50-",thisname, postfix, sep="")
  write.table(thisres, 
            file=difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  print( paste("Diff genes with fdr<0.05 for ", thislabel, " vs DCIS: ",
               sum(thisres$adj.P.Val < 0.05), sep=""))
  
  signchangetab[rownames(thisres)[thisres$adj.P.Val<0.05], thislabel] = "UP"
  signchangetab[rownames(thisres)[thisres$adj.P.Val<0.05], thislabel][thisres$logFC[thisres$adj.P.Val<0.05]<0] = "DOWN" 
}
```

```
## [1] "Diff genes with fdr<0.05 for Normallike vs DCIS: 39"
## [1] "Diff genes with fdr<0.05 for LumA vs DCIS: 31"
## [1] "Diff genes with fdr<0.05 for LumB vs DCIS: 37"
## [1] "Diff genes with fdr<0.05 for Her2 vs DCIS: 17"
## [1] "Diff genes with fdr<0.05 for Basallike vs DCIS: 58"
```

```r
write.table(signchangetab, file=paste(outputdir, "/signchangetab-pam50subgroups", 
                     "-vs-DCIS", postfix, sep=""),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA, na="")

correlationmatrix1 = matrix(nrow=length(pam50labels), ncol=length(pam50labels),dimnames=list(pam50labels,pam50labels))
correlationmatrix2 = correlationmatrix1
for(i in 1:(length(pam50labels)-1))
{
  for(j in (i+1):length(pam50labels))
  {
    correlationmatrix1[i,j] = sum(signchangetab[, pam50labels[i]]==signchangetab[, pam50labels[j]], na.rm=TRUE)
    correlationmatrix1[j,i] = correlationmatrix1[i,j]
    
    correlationmatrix2[i,j] = sum(signchangetab[, pam50labels[i]]!=signchangetab[, pam50labels[j]], na.rm=TRUE)
    correlationmatrix2[j,i] = correlationmatrix2[i,j]
  }

}

fn = paste(outputdir, "/subtypecorrelation", postfix, sep="")
write("Count of significant genes found going in the same direction for the pam50 subtypes", file=fn)
write.table(correlationmatrix1, file=fn,quote=FALSE, sep="\t", col.names=NA, na="X", append=TRUE)
```

```
## Warning in write.table(correlationmatrix1, file = fn, quote = FALSE, sep =
## "\t", : appending column names to file
```

```r
write("\n\n\nCount of significant genes found going in the opposite direction for the pam50 subtypes", file=fn, append=TRUE)
write.table(correlationmatrix2, file=fn,quote=FALSE, sep="\t", col.names=NA,  na="X", append=TRUE)
```

```
## Warning in write.table(correlationmatrix2, file = fn, quote = FALSE, sep =
## "\t", : appending column names to file
```

This is diff genes between all the DCIS (which is a small and diverse group) and each of the pam50 classes for the invasive samples.

<br/>
<br/>


## Stratification based on Immunohistochemistry (IHC) classification.

Many of the samples are also diagnosed based on immunohistochemistry of HER2, PGR and ER proteins. Four diagnoses are used, 
"HER2neg/ERneg/PGRneg", "HER2neg/ERpos", "HER2pos/ERneg", "HER2pos/ERpos". 
We would like too look for microRNA that are changed between all the DCIS and the different IHC diagnoses of invasive samples.


```r
sampleannotation$IHC[sampleannotation$IHC==""] = "Unknown"
sampleannotation$IHC[sampleannotation$tissue_type=="DCIS"] = "DCIS"
table(sampleannotation[,c("tissue_type", "IHC")], useNA="ifany")
```

```
##            IHC
## tissue_type DCIS HER2neg_ERneg_PGRneg HER2neg_ERpos HER2pos_ERneg
##    benign      0                    5             2             0
##    DCIS       18                    0             0             0
##    invasive    0                  199           940            87
##    normal      0                    0             0             0
##            IHC
## tissue_type HER2pos_ERpos Unknown
##    benign               0      24
##    DCIS                 0       0
##    invasive            93      19
##    normal               0     186
```


The diff test,

```r
IHClabels = unique(sampleannotation$IHC[!sampleannotation$IHC %in% c("", "Unknown", "DCIS")])
IHCcontrasts = paste(IHClabels, "-DCIS", sep="")


group = factor(sampleannotation$IHC[dcis_invasive])
batch = factor(sampleannotation$provider[dcis_invasive])

#simple permutation sanity check
# a = batch=="AHUS" & group %in% c("DCIS", IHClabels)
# group[a] = sample(group[a])
# b = batch=="UCAM" & group %in% c("DCIS", IHClabels)
# group[b] = sample(group[b])

design = model.matrix(~0+group+batch)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix[,dcis_invasive], design)
cont.matrix = makeContrasts ( contrasts=IHCcontrasts, levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

for(thiscontrast in IHCcontrasts)
{
  thisres = topTable(fit2, coef=thiscontrast, adjust="BH", number=999)

  thisname= gsub("-", "-vs-", thiscontrast)
  thislabel = strsplit(thiscontrast, "-")[[1]][1]

  difflistfile=paste(outputdir, "/difflist-IHC-",thisname, postfix, sep="")
  write.table(thisres, 
            file=difflistfile,
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  print( paste("Diff genes with fdr<0.05 for ", thislabel, " vs DCIS: ",
               sum(thisres$adj.P.Val < 0.05), sep=""))
  
}
```

```
## [1] "Diff genes with fdr<0.05 for HER2neg_ERpos vs DCIS: 13"
## [1] "Diff genes with fdr<0.05 for HER2pos_ERneg vs DCIS: 21"
## [1] "Diff genes with fdr<0.05 for HER2pos_ERpos vs DCIS: 9"
## [1] "Diff genes with fdr<0.05 for HER2neg_ERneg_PGRneg vs DCIS: 31"
```



**Important note regarding Limma results**
For the diff-tests I ran a few sanity checks where I permuted the group labels whitin batches. For every 3 of 10 run I got a few genes, and once in a while I did get many more. This indicates that some assumtions of the model might be broken to a degree that the resulting p-values should be considered somewhat optimistic.



## References

Breast cancer signatures for invasiveness and prognosis defined by deep sequencing of microRNA.  
Volinia S, Galasso M, Sana ME, Wise TF, Palatini J, Huebner K, Croce CM.  
Proc Natl Acad Sci U S A. 2012 Feb 21;109(8):3024-9. doi: 10.1073/pnas.1200010109. Epub 2012 Feb 6.

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
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] sva_3.10.0            mgcv_1.8-3            nlme_3.1-118         
 [4] corpcor_1.6.7         RColorBrewer_1.0-5    knitr_1.8            
 [7] AgiMicroRna_2.14.0    affycoretools_1.36.1  GO.db_2.14.0         
[10] RSQLite_1.0.0         DBI_0.3.1             AnnotationDbi_1.26.1 
[13] GenomeInfoDb_1.0.2    preprocessCore_1.26.1 affy_1.42.3          
[16] limma_3.20.9          Biobase_2.24.0        BiocGenerics_0.10.0  
[19] xtable_1.7-4         

loaded via a namespace (and not attached):
 [1] acepack_1.3-3.3           affyio_1.32.0            
 [3] annaffy_1.36.0            annotate_1.42.1          
 [5] AnnotationForge_1.6.1     base64enc_0.1-2          
 [7] BatchJobs_1.4             BBmisc_1.7               
 [9] BiocInstaller_1.14.3      BiocParallel_0.6.1       
[11] biomaRt_2.20.0            Biostrings_2.32.1        
[13] biovizBase_1.12.3         bit_1.1-12               
[15] bitops_1.0-6              brew_1.0-6               
[17] BSgenome_1.32.0           Category_2.30.0          
[19] caTools_1.17.1            checkmate_1.5.0          
[21] cluster_1.15.3            codetools_0.2-9          
[23] colorspace_1.2-4          DESeq2_1.4.5             
[25] dichromat_2.0-0           digest_0.6.4             
[27] edgeR_3.6.8               evaluate_0.5.5           
[29] fail_1.2                  ff_2.2-13                
[31] foreach_1.4.2             foreign_0.8-61           
[33] formatR_1.0               Formula_1.1-2            
[35] gcrma_2.36.0              gdata_2.13.3             
[37] genefilter_1.46.1         geneplotter_1.42.0       
[39] GenomicAlignments_1.0.6   GenomicFeatures_1.16.3   
[41] GenomicRanges_1.16.4      ggbio_1.12.10            
[43] ggplot2_1.0.0             GOstats_2.30.0           
[45] gplots_2.14.2             graph_1.42.0             
[47] grid_3.1.1                gridExtra_0.9.1          
[49] GSEABase_1.26.0           gtable_0.1.2             
[51] gtools_3.4.1              Hmisc_3.14-5             
[53] htmltools_0.2.6           hwriter_1.3.2            
[55] IRanges_1.22.10           iterators_1.0.7          
[57] KernSmooth_2.23-13        lattice_0.20-29          
[59] latticeExtra_0.6-26       locfit_1.5-9.1           
[61] markdown_0.7.4            MASS_7.3-35              
[63] Matrix_1.1-4              mime_0.2                 
[65] munsell_0.4.2             nnet_7.3-8               
[67] oligoClasses_1.26.0       PFAM.db_2.14.0           
[69] plyr_1.8.1                proto_0.3-10             
[71] R.methodsS3_1.6.1         R.oo_1.18.0              
[73] R.utils_1.34.0            R2HTML_2.3.1             
[75] RBGL_1.40.1               Rcpp_0.11.3              
[77] RcppArmadillo_0.4.450.1.0 RCurl_1.95-4.3           
[79] ReportingTools_2.4.0      reshape2_1.4             
[81] rmarkdown_0.3.3           rpart_4.1-8              
[83] Rsamtools_1.16.1          rtracklayer_1.24.2       
[85] scales_0.2.4              sendmailR_1.2-1          
[87] splines_3.1.1             stats4_3.1.1             
[89] stringr_0.6.2             survival_2.37-7          
[91] tools_3.1.1               VariantAnnotation_1.10.5 
[93] XML_3.98-1.1              XVector_0.4.0            
[95] zlibbioc_1.10.0          
```

generation ended 2014-12-10 16:24:47. Time spent 0 minutes .

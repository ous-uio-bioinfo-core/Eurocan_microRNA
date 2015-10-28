Consensus, validation, curation.
========================================================
2015-10-28 12:45:21


<br/>
<br/>

We will check for

- concordance between our two approaches called **meta** and **merged**
- concordance between our results and Volinia et al.
- microRNA evidence status for the results

And make most of the tables and plots used in our article.

<BR/>
<BR/>
Read in prepared data.


```r
library(xtable)
if(!exists("inputisread"))
	source("read_input.r")
postfix = ".csv"
separator=","
options(digits=4)
outputdir = "output-analysis"
divdir = paste(outputdir, "/div", sep="")
if(!file.exists(divdir))
	dir.create(divdir)
articledir = paste(outputdir, "/used_in_article", sep="")
if(!file.exists(articledir))
	dir.create(articledir)
supplementaryfile1dir = paste(articledir, "/additionalfile1", sep="")
if(!file.exists(supplementaryfile1dir))
	dir.create(supplementaryfile1dir)
consensusresultsdir=paste(supplementaryfile1dir, "/consensuslists", sep="")
if(!file.exists(consensusresultsdir))
	dir.create(consensusresultsdir)
iClustdifflistdir = paste(supplementaryfile1dir, "/difflists_iClust_UCAM_only", sep="")
if(!file.exists(iClustdifflistdir))
	dir.create(iClustdifflistdir)

fdrCO = 0.05
```
<br/>
<br/>

## Meta vs Merged results.

Our data sets came from two providers, "**AHUS**" and "**METABRIC**" (sometimes called "UCAM"), using two different Agilent microRNA platforms. Batch effects were present as visualized in the QC report. We tested for differentially expressed genes using two approaches, combining all data in Limma using provider as a blocking factor ("**merged**" approach), and  analyzing the two sets independently with Limma and combine based on the p-values ("**meta**" approach).

Here we will list the agreement between the two approaches.

First a recap of the number of samples and grouping based on annotation.

Sample count based on tissue type and provider.


```r
print(xtable(table(sampleannotation[, c("provider" ,"tissue_type")]), 
             caption="", digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> benign </th> <th> DCIS </th> <th> invasive </th> <th> normal </th>  </tr>
  <tr> <td align="right"> AHUS </td> <td align="right">   23 </td> <td align="right">    8 </td> <td align="right">   55 </td> <td align="right">   70 </td> </tr>
  <tr> <td align="right"> UCAM </td> <td align="right">    8 </td> <td align="right">   10 </td> <td align="right"> 1283 </td> <td align="right">  116 </td> </tr>
   </table>
<br/>
<br/>
Sample count based on PAM50 classification in each tissue and from each provider:


```r
for(p in unique(sampleannotation$provider))
{

	print(xtable(table(sampleannotation[sampleannotation$provider==p &
																			sampleannotation$tissue_type %in% c("DCIS", "invasive")	, 
																			c("tissue_type" ,"pam50_est")]), 
             caption=p, digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
}
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> Basallike </th> <th> DCIS </th> <th> Her2 </th> <th> LumA </th> <th> LumB </th> <th> Normallike </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    5 </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">   16 </td> <td align="right">   14 </td> <td align="right">    8 </td> <td align="right">    4 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> UCAM </caption>
<tr> <th>  </th> <th> Basallike </th> <th> DCIS </th> <th> Her2 </th> <th> LumA </th> <th> LumB </th> <th> Normallike </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">   10 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">  193 </td> <td align="right">    0 </td> <td align="right">  159 </td> <td align="right">  460 </td> <td align="right">  370 </td> <td align="right">   90 </td> <td align="right">   11 </td> </tr>
   </table>

Due to the low number of DCIS samples, the pam50 classification were disregarded for those. The comparisons made, were between DCIS samples and the different pam50 classifications of invasive samples. Benign and normal samples were not used in these comparisons.
<br/>
<br/>
Sample count based on IHC diagnosis in each tissue and from each provider


```r
for(p in unique(sampleannotation$provider))
{

	print(xtable(table(sampleannotation[sampleannotation$provider==p , c("tissue_type" ,"IHC")]), 
             caption=p, digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
}
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> HER2neg_ERneg_PGRneg </th> <th> HER2neg_ERpos </th> <th> HER2pos_ERneg </th> <th> HER2pos_ERpos </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> benign </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   23 </td> </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    8 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    1 </td> <td align="right">   28 </td> <td align="right">    4 </td> <td align="right">    8 </td> <td align="right">   14 </td> </tr>
  <tr> <td align="right"> normal </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   70 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> UCAM </caption>
<tr> <th>  </th> <th> HER2neg_ERneg_PGRneg </th> <th> HER2neg_ERpos </th> <th> HER2pos_ERneg </th> <th> HER2pos_ERpos </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> benign </td> <td align="right">    5 </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    2 </td> <td align="right">    3 </td> <td align="right">    5 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">  198 </td> <td align="right">  912 </td> <td align="right">   83 </td> <td align="right">   85 </td> <td align="right">    5 </td> </tr>
  <tr> <td align="right"> normal </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">  116 </td> </tr>
   </table>

Again, due to the low number of DCIS samples, the IHC diagnoses were disregarded for those. The comparisons made, were between DCIS samples and the different IHC diagnoses of invasive samples. Benign and normal samples were not used in these comparisons either.


<br/>
<br/>

```r
# creating table1 our in article.
notbenign = !sampleannotation$tissue_type=="benign"
invasiveonly = sampleannotation$tissue_type=="invasive"
a=table(sampleannotation[notbenign, c("tissue_type", "provider")])
b=table(sampleannotation[invasiveonly, c("pam50_est", "provider")])
c=table(sampleannotation[invasiveonly, c("IHC", "provider")])
table1 = rbind(a[c("normal", "DCIS", "invasive"),], 
							 b[c("LumA", "LumB", "Her2", "Basallike", "Normallike"),],
							 c[c("HER2neg_ERpos", "HER2pos_ERpos", "HER2pos_ERneg", "HER2neg_ERneg_PGRneg"),],
							 total=colSums(a))
write.table(table1, file=paste(articledir, "/table1", postfix, sep=""), sep=separator, col.names=NA, quote=FALSE)
```

Read in our results from files.


```r
#metafilefolder = "not_in_github/lilianas"
metafilefolder = "output-analysis/difflists_meta"
mergedfilefolder = "output-analysis/difflists_merged"

meta2merged = data.frame(metafn=character(0), mergedfn=character(0), direction=numeric(0), comparisontype=character(0), stringsAsFactors=FALSE)
meta2merged["DCIS-normal",] = list("DCIS_to_normal_limma.txt", "difflist-DCIS-vs-normal.txt",-1, "tissue")
meta2merged["invasive-DCIS",] = list("DCIS_to_invasive_limma.txt", "difflist-invasive-vs-DCIS.txt",1, "tissue")
meta2merged["Normallike-DCIS",] = list("dcis_to_normallike_limma.txt", "difflist-pam50-Normallike-vs-DCIS.txt",1, "pam50")
meta2merged["LumA-DCIS",] = list("dcis_to_lumA_limma.txt", "difflist-pam50-LumA-vs-DCIS.txt",1, "pam50")
meta2merged["LumB-DCIS",] = list("dcis_to_lumB_limma.txt", "difflist-pam50-LumB-vs-DCIS.txt",1, "pam50")
meta2merged["Her2-DCIS",] = list("dcis_to_her2_limma.txt", "difflist-pam50-Her2-vs-DCIS.txt",1, "pam50")
meta2merged["Basallike-DCIS",] = list("basallike_to_dcis_limma.txt", "difflist-pam50-Basallike-vs-DCIS.txt",-1, "pam50")
meta2merged["HER2neg_ERneg_PGRneg-DCIS",] = list("dcis_to_HER2neg_ERneg_PGRneg_limma.txt", "difflist-IHC-HER2neg_ERneg_PGRneg-vs-DCIS.txt",1, "IHC")
meta2merged["HER2neg_ERpos-DCIS",] = list("dcis_to_HER2neg_ERpos_limma.txt", "difflist-IHC-HER2neg_ERpos-vs-DCIS.txt",1, "IHC")
meta2merged["HER2pos_ERneg-DCIS",] = list("dcis_to_HER2pos_ERneg_limma.txt", "difflist-IHC-HER2pos_ERneg-vs-DCIS.txt",1, "IHC")
meta2merged["HER2pos_ERpos-DCIS",] = list("dcis_to_HER2pos_ERpos_limma.txt", "difflist-IHC-HER2pos_ERpos-vs-DCIS.txt",1, "IHC")

metatables = list()
mergedtables = list()
consensustables = list()
summarytab = data.frame( found_meta=integer(0),  found_merged=integer(0), 
                       found_both=integer(0),  stringsAsFactors=FALSE)

for(i in 1:nrow(meta2merged))
{   
    thiscomp = rownames(meta2merged)[i]
    fn = paste(metafilefolder, "/", meta2merged[thiscomp, "metafn"], sep="")
    if(file.info(fn)$size>10)
    	{
    thismetatab = read.table(fn, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="")
    } else{
    	thismetatab = thismetatab[0,] # workaround to deal with empty files from meta.
    }
    fn = paste(mergedfilefolder, "/", meta2merged[thiscomp, "mergedfn"], sep="")
    thismergedtab = read.table(fn, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="", row.names=1)
    
    #harmonizing a little
    rownames(thismetatab) = thismetatab$MIMAT
    thismetatab$TestStatistic = thismetatab$TestStatistic * meta2merged[thiscomp, "direction"]
  #summarytab[thiscomp, "metafile"] = rownames(thisrows)[1]
  #summarytab[thiscomp, "mergedfile"] = thisrows[1, "mergedfn"]

  a = thismetatab$adj.P.Val <= 0.05
  summarytab[thiscomp, "found_meta"] =sum(a)
  b = thismergedtab$adj.P.Val <= 0.05
  summarytab[thiscomp, "found_merged"] = sum(b)
  
  commonMIMATID = intersect( rownames(thismetatab)[a],  rownames(thismergedtab)[b])
  summarytab[thiscomp, "found_both"] = sum( thismetatab[commonMIMATID, "TestStatistic"] * thismergedtab[commonMIMATID, "logFC"] > 0) 
  summarytab[thiscomp, "wrong_direction"] = sum( thismetatab[commonMIMATID, "TestStatistic"] * thismergedtab[commonMIMATID, "logFC"] < 0) 
  rm(commonMIMATID)
  metatables[[thiscomp]] = thismetatab
  mergedtables[[thiscomp]] = thismergedtab
  
  x = mergedtables[[thiscomp]][, c("logFC", "adj.P.Val")]
  y = metatables[[thiscomp]][, c("TestStatistic", "adj.P.Val")]
  names(x) = paste(names(x), ".merged", sep="")
  names(y) = paste(names(y), ".meta", sep="")
  
 	thistable = merge(x, y, by.x = "row.names", by.y="row.names",  all.x=TRUE)
  row.names(thistable) = thistable$Row.names
  thistable = cbind(MIMAT = thistable$Row.names, name=miRNA2MIMAT[thistable$Row.names, "preferredname"], thistable[,-1])
  thistable$same_dir =  (thistable$logFC.merged * thistable$TestStatistic.meta ) > 0
  thistable$consensus = thistable$adj.P.Val.merged <= 0.05 & 
  	thistable$adj.P.Val.meta <= 0.05 & thistable$same_dir
  consensustables[[thiscomp]] = thistable[order(thistable$consensus, -thistable$adj.P.Val.merged, decreasing=TRUE),]
}

if(sum(summarytab$wrong_direction)>0)
	warning("Some common microRNAs were found that did not have the same dirrection of change!!")
```
The above tables are written to [files](output-analysis/used_in_article/additionalfile1/consensuslists) later.
<br/>
<br/>
Next, summarize the number of microRNAs found differentially expressed between states, according to the two approaches. The direction of fold change is taken into account. 


```r
	print(xtable(summarytab[,1:3], 
             caption="Overlap between meta and merged", digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:22 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> Overlap between meta and merged </caption>
<tr> <th>  </th> <th> found_meta </th> <th> found_merged </th> <th> found_both </th>  </tr>
  <tr> <td align="right"> DCIS-normal </td> <td align="right">  116 </td> <td align="right">   86 </td> <td align="right">   82 </td> </tr>
  <tr> <td align="right"> invasive-DCIS </td> <td align="right">    0 </td> <td align="right">   10 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Normallike-DCIS </td> <td align="right">    9 </td> <td align="right">   37 </td> <td align="right">    9 </td> </tr>
  <tr> <td align="right"> LumA-DCIS </td> <td align="right">    9 </td> <td align="right">   30 </td> <td align="right">    6 </td> </tr>
  <tr> <td align="right"> LumB-DCIS </td> <td align="right">   30 </td> <td align="right">   44 </td> <td align="right">   24 </td> </tr>
  <tr> <td align="right"> Her2-DCIS </td> <td align="right">   11 </td> <td align="right">   16 </td> <td align="right">    6 </td> </tr>
  <tr> <td align="right"> Basallike-DCIS </td> <td align="right">    7 </td> <td align="right">   62 </td> <td align="right">    7 </td> </tr>
  <tr> <td align="right"> HER2neg_ERneg_PGRneg-DCIS </td> <td align="right">    2 </td> <td align="right">   30 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> HER2neg_ERpos-DCIS </td> <td align="right">    1 </td> <td align="right">   11 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERneg-DCIS </td> <td align="right">   15 </td> <td align="right">   19 </td> <td align="right">    8 </td> </tr>
  <tr> <td align="right"> HER2pos_ERpos-DCIS </td> <td align="right">    3 </td> <td align="right">    5 </td> <td align="right">    0 </td> </tr>
   </table>

```r
#	write.table(summarytab[,1:3], file=paste(supplementaryfile1dir, "/suppltable1", postfix, sep=""), 
#							row.names=TRUE, col.names=NA, sep=separator)
```


For several of the comparisons, the shorter list seems to be mostly included in the longer. We will use the intersect as differentially expressed microRNAs, i.e a microRNA is deemed differentially expressed if it is within 0.05 FDR in both the meta and merged analyses.
<br/>
<br/>

Next, a table is written to file with all the expressed (not filtered out) microRNAs listed and whether they were found to be up regulated or down regulated in the different comparisons.


```r
signchangetab = as.data.frame(matrix( nrow=nrow(consensustables[[1]]), 
																			ncol=length(consensustables)+2), stringsAsFactors=FALSE)
colnames(signchangetab) = c("MIMATID", "name", names(consensustables))
rownames(signchangetab) = rownames(consensustables[[1]])
signchangetab[,"MIMATID"] = consensustables[[1]]$MIMAT
signchangetab[,"name"] = miRNA2MIMAT[consensustables[[1]]$MIMAT, "MIMAT"]

for(n in names(consensustables))
{
	x = consensustables[[n]][rownames(signchangetab),]
	signchangetab[ x$consensus & x$logFC.merged<0, n ] = "DOWN"
	signchangetab[ x$consensus & x$logFC.merged>0, n ] = "UP" 
}
signchangetabfn = paste(divdir, "/signchangetab", postfix, sep="")
write.table(signchangetab, file=signchangetabfn,
            quote=FALSE, sep=separator, row.names=FALSE, col.names=TRUE, na="")
```

[Link to file](output-analysis/div/signchangetab.csv)
<br/>
<br/>
Next, the overlap in findings between the comparisons.


```r
comps = names(consensustables)
correlationmatrix1 = matrix(nrow=length(comps), ncol=length(comps),dimnames=list(comps,comps))
correlationmatrix2 = correlationmatrix1
for(i in 1:(length(comps)-1))
{
  for(j in (i+1):length(comps))
  {
    correlationmatrix1[i,j] = sum(signchangetab[, comps[i]]==signchangetab[, comps[j]], na.rm=TRUE)
    correlationmatrix1[j,i] = correlationmatrix1[i,j]
    
    correlationmatrix2[i,j] = sum(signchangetab[, comps[i]]!=signchangetab[, comps[j]], na.rm=TRUE)
    correlationmatrix2[j,i] = correlationmatrix2[i,j]
  }

}

correlationmatrixfn = paste(divdir, "/subtypecorrelation", postfix, sep="")
h1 = "Count of significant microRNAs found going in the same direction for the subtypes"
write(h1, file=correlationmatrixfn)
write(paste(separator, colnames(correlationmatrix1), collapse=separator), file=correlationmatrixfn, append=TRUE)
write.table(correlationmatrix1, file=correlationmatrixfn,quote=FALSE, sep=separator, col.names=FALSE, na="X", append=TRUE)
write("\n\n\n", file=correlationmatrixfn, append=TRUE)
h2 = "Count of significant microRNAs found going in the opposite direction for the subtypes"
write(h2, file=correlationmatrixfn, append=TRUE)
write(paste(separator, colnames(correlationmatrix2), collapse=separator), file=correlationmatrixfn, append=TRUE)
write.table(correlationmatrix2, file=correlationmatrixfn,quote=FALSE, sep=separator, col.names=FALSE,  na="X", append=TRUE)

print(xtable(correlationmatrix1, 
             caption=h1, digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Oct 28 12:45:23 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> Count of significant microRNAs found going in the same direction for the subtypes </caption>
<tr> <th>  </th> <th> DCIS-normal </th> <th> invasive-DCIS </th> <th> Normallike-DCIS </th> <th> LumA-DCIS </th> <th> LumB-DCIS </th> <th> Her2-DCIS </th> <th> Basallike-DCIS </th> <th> HER2neg_ERneg_PGRneg-DCIS </th> <th> HER2neg_ERpos-DCIS </th> <th> HER2pos_ERneg-DCIS </th> <th> HER2pos_ERpos-DCIS </th>  </tr>
  <tr> <td align="right"> DCIS-normal </td> <td align="right">  </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">   10 </td> <td align="right">    2 </td> <td align="right">    2 </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive-DCIS </td> <td align="right">    0 </td> <td align="right">  </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Normallike-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">  </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> LumA-DCIS </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    2 </td> <td align="right">  </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> LumB-DCIS </td> <td align="right">   10 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    2 </td> <td align="right">  </td> <td align="right">    3 </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    2 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Her2-DCIS </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    3 </td> <td align="right">  </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    4 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Basallike-DCIS </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">  </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2neg_ERneg_PGRneg-DCIS </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">    2 </td> <td align="right">  </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2neg_ERpos-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">  </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERneg-DCIS </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    2 </td> <td align="right">    4 </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">    0 </td> <td align="right">  </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERpos-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">  </td> </tr>
   </table>

[Link to file](output-analysis/div/subtypecorrelation.csv)

<br/>
<br/>
<br/>

## Our results vs. Volinia et al.´s results


```r
volinia_DCIS_vs_normal_file = "volinia-DCIS-vs-normal.txt"
volinia_invasive_vs_DCIS_file = "volinia-invasive-vs-DCIS.txt"
```

volinia-DCIS-vs-normal.txt and  volinia-invasive-vs-DCIS.txt are copy-pasted from Volinia et al.s [supplementary tables S1 and S2](http://www.pnas.org/content/suppl/2012/02/02/1200010109.DCSupplemental/pnas.201200010SI.pdf). The tables were reformatted from the pdf file into a more practical tabular text format. This is the results we aim to validate with the AHUS and METABRIC data sets.

Read in Volinia results, and do a minor re-formatting.


```r
voliniatab = list()
voliniatab[["invasive-DCIS"]] = read.table(file=paste(annotdir, "/", volinia_invasive_vs_DCIS_file, sep=""),
                                      sep="\t", check.names=TRUE, 
                                      stringsAsFactors=FALSE, header=TRUE)
voliniatab[["DCIS-normal"]]  = read.table(file=paste(annotdir, "/", volinia_DCIS_vs_normal_file, sep=""),
                                    sep="\t", check.names=TRUE, 
                                    stringsAsFactors=FALSE, header=TRUE)

for(n in names(voliniatab))
{
	# add a MIMATID 
	voliniatab[[n]] = as.data.frame(append(voliniatab[[n]], 
                  list(MIMATID=miRNA2MIMAT$MIMAT[match(voliniatab[[n]]$miRNA, miRNA2MIMAT$volinia)]), 
                  after=1))
	# use log2 for foldchange in this document.
	voliniatab[[n]]$fc=log2(voliniatab[[n]]$fc)
	names(voliniatab[[n]])[names(voliniatab[[n]])=="fc"] = "log2FC"
	names(voliniatab[[n]])[names(voliniatab[[n]])=="pval"] = "FDR" # this was described as FDR in their paper.
}
```

This is how the tables look now (top only).


```r
print(head(voliniatab[[1]]), row.names = FALSE)
```

```
##       miRNA      MIMATID  log2FC      FDR
##     miR-10b MIMAT0000254 -1.2863 4.00e-07
##     miR-143 MIMAT0000435 -1.3585 5.00e-07
##      let-7d MIMAT0000065  1.3334 5.20e-06
##  miR-218(2) MIMAT0000275 -0.6897 1.59e-05
##  miR-335-5p MIMAT0000765 -1.1844 5.43e-05
##     miR-126 MIMAT0000445 -1.1203 7.13e-05
```


Compared to the original tables from Volinia, the fold change is here log2 transformed and "IDC" is renamed to "invasive", and the microRNAs are also mapped to MIMATID, which will be used in order to map to our results. 


Count the overlaps of results and make a "validated"-call for each microRNA listed by Volinia et al.


```r
volinia_overlap = list()

for(contrast in names(voliniatab))
{
	overlaptab = voliniatab[[contrast]]
	overlaptab = overlaptab[!is.na(overlaptab$MIMATID),]
	rownames(overlaptab) = overlaptab$MIMATID
	overlaptab$Direction = ifelse(overlaptab$log2FC > 0, "UP", "DOWN")
	overlaptab = overlaptab[, c("miRNA", "MIMATID", "Direction","FDR")]
	names(overlaptab) = paste(names(overlaptab), "_vol", sep="")
	
	tab = consensustables[[contrast]]
	a = intersect(  rownames(tab), rownames(overlaptab) )
	overlaptab$Direction_merged = NA
	overlaptab[a, "Direction_merged"] = ifelse(tab[a, "logFC.merged"] > 0,  "UP", "DOWN")
	overlaptab[a, "adj.P.Val.merged"] = tab[a, "adj.P.Val.merged"]
	overlaptab$Direction_meta = NA
	overlaptab[a, "Direction_meta"] = tab[a, "TestStatistic.meta"]
	overlaptab[is.na(overlaptab$Direction_meta), "Direction_meta"] =  0	
	overlaptab[overlaptab$Direction_meta>0, "Direction_meta"] =  "UP"
	overlaptab[overlaptab$Direction_meta<0, "Direction_meta"] =  "DOWN"
	overlaptab[overlaptab$Direction_meta==0, "Direction_meta"] =  NA	
	overlaptab[a, "adj.P.Val.meta"] = tab[a, "adj.P.Val.meta"]
	
	x = rowSums (overlaptab[,grepl("Direction", colnames(overlaptab))] == "UP")
	overlaptab$validated = NA
	overlaptab$validated[x %in% c(0,3)] = TRUE
	overlaptab$validated[x %in% c(1,2)] = FALSE
	overlaptab$validated = overlaptab$validated & (overlaptab$adj.P.Val.merged<=fdrCO & overlaptab$adj.P.Val.meta<=fdrCO)

	overlaptab = overlaptab[ order(overlaptab$adj.P.Val.merged) , ]
	volinia_overlap[[contrast]] = overlaptab
	write.table(overlaptab, 
            file=paste(outputdir, "/Volinia_validation_", contrast, postfix, sep=""),
            quote=FALSE, sep=separator, row.names=TRUE, col.names=NA)

	rm(x, overlaptab, contrast)
}

# write table2 in our article
table2 = consensustables[["DCIS-normal"]]
table2 = table2[table2$consensus,]
table2 = table2[table2$MIMAT %in% rownames( volinia_overlap[["DCIS-normal"]])[volinia_overlap[["DCIS-normal"]]$validated], ]
table2 = table2[order(table2$logFC.merged),2:6]
table2[,c(2,4)] = round(table2[,c(2,4)], 2)
write.table(table2, file=paste(articledir, "/table2", postfix, sep=""), sep=separator, col.names=NA, quote=FALSE)
```

And below is the result list from Volinia et al. with our results attached. First the "DCIS-normal" comparison. Rows are sorted by FDR from our "merged" approach. 


```r
print( xtable( volinia_overlap[["DCIS-normal"]] , caption="microRNA found in Volinia compared to validation data", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> microRNA found in Volinia compared to validation data </caption>
<tr> <th> miRNA_vol </th> <th> MIMATID_vol </th> <th> Direction_vol </th> <th> FDR_vol </th> <th> Direction_merged </th> <th> adj.P.Val.merged </th> <th> Direction_meta </th> <th> adj.P.Val.meta </th> <th> validated </th>  </tr>
  <tr> <td> miR-21 </td> <td> MIMAT0000076 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-497 </td> <td> MIMAT0002820 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-145 </td> <td> MIMAT0000437 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-193a-5p </td> <td> MIMAT0004614 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-378 </td> <td> MIMAT0000732 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-155 </td> <td> MIMAT0000646 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-200c </td> <td> MIMAT0000617 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-425 </td> <td> MIMAT0003393 </td> <td> UP </td> <td align="right"> 0.012 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-183 </td> <td> MIMAT0000261 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-200b </td> <td> MIMAT0000318 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-429 </td> <td> MIMAT0001536 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-96 </td> <td> MIMAT0000095 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-142-3p </td> <td> MIMAT0000434 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-143* </td> <td> MIMAT0004599 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-140-3p </td> <td> MIMAT0004597 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-142-5p </td> <td> MIMAT0000433 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-99a </td> <td> MIMAT0000097 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> let-7b </td> <td> MIMAT0000063 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> let-7c </td> <td> MIMAT0000064 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-193b </td> <td> MIMAT0002819 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-19b(2) </td> <td> MIMAT0000074 </td> <td> DOWN </td> <td align="right"> 0.015 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-106b </td> <td> MIMAT0000680 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-145* </td> <td> MIMAT0004601 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-342 </td> <td> MIMAT0000753 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-100 </td> <td> MIMAT0000098 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.008 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> DOWN </td> <td align="right"> 0.009 </td> <td> UP </td> <td align="right"> 0.022 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-22 </td> <td> MIMAT0000077 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.024 </td> <td> DOWN </td> <td align="right"> 0.049 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-652 </td> <td> MIMAT0003322 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-125b(2) </td> <td> MIMAT0000423 </td> <td> DOWN </td> <td align="right"> 0.004 </td> <td> DOWN </td> <td align="right"> 0.044 </td> <td> DOWN </td> <td align="right"> 0.022 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-374b </td> <td> MIMAT0004955 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.051 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-92a(2) </td> <td> MIMAT0000092 </td> <td> UP </td> <td align="right"> 0.011 </td> <td> DOWN </td> <td align="right"> 0.053 </td> <td> DOWN </td> <td align="right"> 0.020 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-19a </td> <td> MIMAT0000073 </td> <td> UP </td> <td align="right"> 0.008 </td> <td> UP </td> <td align="right"> 0.072 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-376c </td> <td> MIMAT0000720 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> DOWN </td> <td align="right"> 0.088 </td> <td> DOWN </td> <td align="right"> 0.324 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.112 </td> <td> UP </td> <td align="right"> 0.091 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-361-5p </td> <td> MIMAT0000703 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.114 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-30d </td> <td> MIMAT0000245 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.115 </td> <td> DOWN </td> <td align="right"> 0.034 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-376a-3p(2) </td> <td> MIMAT0000729 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.117 </td> <td> DOWN </td> <td align="right"> 0.319 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-223 </td> <td> MIMAT0000280 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> UP </td> <td align="right"> 0.123 </td> <td> UP </td> <td align="right"> 0.079 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-185 </td> <td> MIMAT0000455 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.129 </td> <td> UP </td> <td align="right"> 0.182 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-374a </td> <td> MIMAT0000727 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.194 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-127-3p </td> <td> MIMAT0000446 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.199 </td> <td> DOWN </td> <td align="right"> 0.084 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-26a(2) </td> <td> MIMAT0000082 </td> <td> UP </td> <td align="right"> 0.013 </td> <td> DOWN </td> <td align="right"> 0.235 </td> <td> DOWN </td> <td align="right"> 0.204 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-16(2) </td> <td> MIMAT0000069 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.282 </td> <td> UP </td> <td align="right"> 0.325 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-29b(2) </td> <td> MIMAT0000100 </td> <td> UP </td> <td align="right"> 0.007 </td> <td> UP </td> <td align="right"> 0.290 </td> <td> UP </td> <td align="right"> 0.066 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-26b </td> <td> MIMAT0000083 </td> <td> UP </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.378 </td> <td> UP </td> <td align="right"> 0.266 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-574-3p </td> <td> MIMAT0003239 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.379 </td> <td> UP </td> <td align="right"> 0.271 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-20a </td> <td> MIMAT0000075 </td> <td> UP </td> <td align="right"> 0.007 </td> <td> UP </td> <td align="right"> 0.447 </td> <td> UP </td> <td align="right"> 0.350 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-107 </td> <td> MIMAT0000104 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.469 </td> <td> UP </td> <td align="right"> 0.584 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-125a </td> <td> MIMAT0000443 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.783 </td> <td> UP </td> <td align="right"> 0.628 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-423-5p </td> <td> MIMAT0004748 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.853 </td> <td> UP </td> <td align="right"> 0.756 </td> <td> FALSE </td> </tr>
  <tr> <td> let-7d </td> <td> MIMAT0000065 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.887 </td> <td> DOWN </td> <td align="right"> 0.732 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15b </td> <td> MIMAT0000417 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.945 </td> <td> UP </td> <td align="right"> 0.641 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-29c </td> <td> MIMAT0000681 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.951 </td> <td> DOWN </td> <td align="right"> 0.640 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15a </td> <td> MIMAT0000068 </td> <td> UP </td> <td align="right"> 0.013 </td> <td> DOWN </td> <td align="right"> 0.972 </td> <td> DOWN </td> <td align="right"> 0.981 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-324 </td> <td> MIMAT0000761 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.973 </td> <td> DOWN </td> <td align="right"> 0.583 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-182 </td> <td> MIMAT0000259 </td> <td> UP </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-423-3p </td> <td> MIMAT0001340 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-106b* </td> <td> MIMAT0004672 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-32 </td> <td> MIMAT0000090 </td> <td> UP </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-28-3p </td> <td> MIMAT0004502 </td> <td> DOWN </td> <td align="right"> 0.004 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-340 </td> <td> MIMAT0004692 </td> <td> UP </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-452 </td> <td> MIMAT0001635 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-17* </td> <td> MIMAT0000071 </td> <td> DOWN </td> <td align="right"> 0.015 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
   </table>

<br/>
<br/>
And here is the same table for the "invasive-DCIS" comparison. Both tables are also written to [files](output-analysis).


```r
x = order(volinia_overlap[["invasive-DCIS"]]$adj.P.Val.merged)
print( xtable( volinia_overlap[["invasive-DCIS"]][x,] , caption="microRNA found in Volinia compared to validation data", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> microRNA found in Volinia compared to validation data </caption>
<tr> <th> miRNA_vol </th> <th> MIMATID_vol </th> <th> Direction_vol </th> <th> FDR_vol </th> <th> Direction_merged </th> <th> adj.P.Val.merged </th> <th> Direction_meta </th> <th> adj.P.Val.meta </th> <th> validated </th>  </tr>
  <tr> <td> let-7d </td> <td> MIMAT0000065 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.057 </td> <td> UP </td> <td align="right"> 0.710 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.299 </td> <td> UP </td> <td align="right"> 0.936 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.306 </td> <td> DOWN </td> <td align="right"> 0.710 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-126 </td> <td> MIMAT0000445 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.359 </td> <td> DOWN </td> <td align="right"> 0.821 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-143 </td> <td> MIMAT0000435 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.497 </td> <td> DOWN </td> <td align="right"> 0.959 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-218(2) </td> <td> MIMAT0000275 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.501 </td> <td> DOWN </td> <td align="right"> 0.936 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-181a(2) </td> <td> MIMAT0000256 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.600 </td> <td> UP </td> <td align="right"> 0.977 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-335-5p </td> <td> MIMAT0000765 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.614 </td> <td> DOWN </td> <td align="right"> 0.959 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-10b </td> <td> MIMAT0000254 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.991 </td> <td> UP </td> <td align="right"> 0.959 </td> <td> FALSE </td> </tr>
   </table>
<br/>
<br/>


<br/>
<br/>

## Checking results towards a curated list of microRNA´s

Bastian Fromm and colleagues have carefully evaluated the evidence for the Human micro RNAs reported in mirBase and provided a classification based on this. In short, a lot of microRNAs in mirBASE and included in the microarray platforms are probably not real microRNAs.

A list of human microRNAs and their curation classification for our microRNAs is provided:


```r
	curatedfn = "curated_microRNA_MIMAT.csv"
  orgcuratedtab  = read.table(paste(annotdir, "/", curatedfn, sep=""), sep="\t", 
  												 header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="")
	orgcuratedtab$miRNA = paste("hsa-", orgcuratedtab$miRNA, sep="")
	
	# rearrange to a more practical object
	curatedtab = orgcuratedtab[,c(1,2,3)]
	names(curatedtab)[3]= "MIMAT"
	
	curatedtab = rbind(curatedtab, 
										 setNames( orgcuratedtab[,c(1,2,4)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,5)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,6)], names(curatedtab)) )
	curatedtab = curatedtab[curatedtab$MIMAT != "", ]

	table(orgcuratedtab$status)
```

```
## 
##   purple rejected    white   yellow 
##      302      969      484        4
```

The above counts are the microRNAs that were mapped to MIMAT IDs. 

The status labels explained:

- *purple*: equivocal - expression data lacking from one arm
- *rejected*: not a microRNA based on structure/read distribution
- *white*: a microRNA
- *yellow*: non-canonical microRNA

We want to check the curation status of our findings. Are our significant microRNAs real microRNAs? In addition as a quality control, it is useful to see the proportion of the statuses in the microRNAs that are filtered out and the ones retained.

Adding the curated status to the results, first the microRNAs found in Volinia et al.


```r
for(n in names(volinia_overlap))
{
	x = unlist(lapply(rownames(volinia_overlap[[n]]), FUN=function(x)strsplit(x, "_")[[1]][1]))
	volinia_overlap[[n]]$curated = curatedtab[ match( x, curatedtab$MIMAT ), "status"]
	print( paste( n , names(table(volinia_overlap[[n]]$curated)), table(volinia_overlap[[n]]$curated)))
}
```

```
## [1] "invasive-DCIS white 9"
## [1] "DCIS-normal white 63"
```
The microRNAs reported by Volinia et al. are all in the "for sure microRNA" list in Fromm et al.
<br/>
<br/>
And now the same assignment and summary for our results from the AHUS and METABRIC data sets.


```r
namelists = list()
for(n in names(consensustables))
{
	x = unlist(lapply(rownames(consensustables[[n]]), FUN=function(x)strsplit(x, "_")[[1]][1])) # not neede anymore?
	consensustables[[n]]$curated = curatedtab[ match( x, curatedtab$MIMAT ), "status"]
	write.table( format( consensustables[[n]], digits=4), 
            file=paste(consensusresultsdir, "/consensusresults_",meta2merged[n, "comparisontype"], "__", n, postfix, sep=""),
            quote=FALSE, sep=separator, row.names=FALSE, col.names=TRUE)
	
	# for table3 and 4.
	dd = consensustables[[n]]
	dd = dd[dd$consensus & !is.na(dd$curated),]
	dd = dd[order(dd$logFC.merged, decreasing=TRUE), ]
	v = as.character(dd$name)
	v[dd$curated!="white"] = paste("(", v[dd$curated!="white"], ")", sep="")
	v[dd$logFC.merged>0] = paste( v[dd$logFC.merged>0], " +", sep="")
	v[dd$logFC.merged<0] = paste( v[dd$logFC.merged<0], " -", sep="")
	namelists[[n]]=v
	
}

# add a nice version of the sampleannotation to the suppl.
write.table(sampleannotation[,c("sample_id", "tissue_type", "pam50_est", "IHC", "iClust", "provider")], 
						file=paste(supplementaryfile1dir, "/sampleannotation_suppl", postfix, sep=""), sep=separator, 
						row.names=FALSE, col.names=TRUE)

# copy files made elsewhere to be included in suppl.
file.copy(from =list.files("input/additionalfile", full.names=TRUE), to=supplementaryfile1dir, overwrite=TRUE)
```

[1] TRUE TRUE TRUE TRUE

```r
# write table3 and table4 in our article
tableclasses = c("LumA-DCIS", "LumB-DCIS", "Her2-DCIS", "Basallike-DCIS", "Normallike-DCIS")
nrows = max( unlist(lapply( namelists[tableclasses], FUN=length) ))
table3 = matrix("",ncol=length(tableclasses), nrow=nrows)
colnames(table3)=tableclasses
for(n in tableclasses)
	table3[1:length(namelists[[n]]), n] = namelists[[n]]
write.table(table3, file=paste(articledir, "/table3", postfix, sep=""), sep=separator, row.names=FALSE, quote=FALSE)

tableclasses = c("HER2pos_ERpos-DCIS", "HER2neg_ERpos-DCIS", "HER2pos_ERneg-DCIS", "HER2neg_ERneg_PGRneg-DCIS")
nrows = max( unlist(lapply( namelists[tableclasses], FUN=length) ))
table4 = matrix("",ncol=length(tableclasses), nrow=nrows)
colnames(table4)=tableclasses
for(n in tableclasses)
	if(length(namelists[[n]])>0)
		table4[1:length(namelists[[n]]), n] = namelists[[n]]
write.table(table4, file=paste(articledir, "/table4", postfix, sep=""), sep=separator, row.names=FALSE, quote=FALSE)
```

Result files for all the comparisons with p-values for the two approaches, and with the curation status added, are [here](output-analysis/used_in_article/additionalfile1/consensuslists). The tables look like this (top of the DCIS-normal comparison as an example):


```r
print( xtable( consensustables[["DCIS-normal"]][1:5,], caption="", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> MIMAT </th> <th> name </th> <th> logFC.merged </th> <th> adj.P.Val.merged </th> <th> TestStatistic.meta </th> <th> adj.P.Val.meta </th> <th> same_dir </th> <th> consensus </th> <th> curated </th>  </tr>
  <tr> <td align="right"> MIMAT0005893 </td> <td> MIMAT0005893 </td> <td> hsa-miR-1305 </td> <td align="right"> 0.681 </td> <td align="right"> 0.000 </td> <td align="right"> 6.120 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> purple </td> </tr>
  <tr> <td align="right"> MIMAT0005942 </td> <td> MIMAT0005942 </td> <td> hsa-miR-1288-3p </td> <td align="right"> 0.342 </td> <td align="right"> 0.000 </td> <td align="right"> 5.095 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> rejected </td> </tr>
  <tr> <td align="right"> MIMAT0000076 </td> <td> MIMAT0000076 </td> <td> hsa-miR-21-5p </td> <td align="right"> 2.407 </td> <td align="right"> 0.000 </td> <td align="right"> 8.676 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td align="right"> MIMAT0000250 </td> <td> MIMAT0000250 </td> <td> hsa-miR-139-5p </td> <td align="right"> -1.276 </td> <td align="right"> 0.000 </td> <td align="right"> -8.167 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td align="right"> MIMAT0005938 </td> <td> MIMAT0005938 </td> <td> hsa-miR-1274b </td> <td align="right"> 1.600 </td> <td align="right"> 0.000 </td> <td align="right"> 6.308 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td>  </td> </tr>
   </table>

**Curation status for the significant microRNAs**
Below is a table with the counts of microRNAs in each curation class for the comparisons.

```r
curclasses = unique(curatedtab$status)
curclasscounts = matrix(nrow=length(consensustables), ncol=length(curclasses)+1)
rownames(curclasscounts) = names(consensustables)
colnames(curclasscounts) = c(curclasses, "NA")
for(n in names(consensustables))
{
	a = consensustables[[n]]$consensus
	for(thisclass in curclasses)
		curclasscounts[n,thisclass] = sum(consensustables[[n]]$curated[a]==thisclass, na.rm =TRUE)
	curclasscounts[n,"NA"] = sum(is.na(consensustables[[n]]$curated[a]))
}
print( xtable( curclasscounts, caption="", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> white </th> <th> yellow </th> <th> purple </th> <th> rejected </th> <th> NA </th>  </tr>
  <tr> <td align="right"> DCIS-normal </td> <td align="right">   68 </td> <td align="right">    1 </td> <td align="right">    1 </td> <td align="right">    9 </td> <td align="right">    3 </td> </tr>
  <tr> <td align="right"> invasive-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Normallike-DCIS </td> <td align="right">    5 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    3 </td> </tr>
  <tr> <td align="right"> LumA-DCIS </td> <td align="right">    6 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> LumB-DCIS </td> <td align="right">   24 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Her2-DCIS </td> <td align="right">    6 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Basallike-DCIS </td> <td align="right">    7 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2neg_ERneg_PGRneg-DCIS </td> <td align="right">    2 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2neg_ERpos-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERneg-DCIS </td> <td align="right">    7 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERpos-DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
   </table>


<br/>
<br/>

**Differences in curation status for groups of microRNA**

The difference in curated status for the filtered out microRNAs and the ones kept as expressed, and the ones found as significantly differentially expressed between the "DCIS" and "normal" state are interesting at least from a quality control point of view.



```r
filteredcommonMIMAT = consensustables[[1]]$MIMAT
unfilteredcommonMIMAT = miRNA2MIMAT$MIMAT[ !is.na(miRNA2MIMAT$UCAM) &
																		 !is.na(miRNA2MIMAT$AHUS) & 
																		 !miRNA2MIMAT$MIMAT %in% filteredcommonMIMAT]

filt1 = unlist(lapply(filteredcommonMIMAT, FUN=function(x)strsplit(x, "_")[[1]][1]))
unfilt1 = unlist(lapply(unfilteredcommonMIMAT, FUN=function(x)strsplit(x, "_")[[1]][1]))

expressed =	table(curatedtab$status[match(filt1, curatedtab$MIMAT)])
nonexpressed=	table(curatedtab$status[match(unfilt1[!unfilt1 %in% filt1],curatedtab$MIMAT) ])

significant = table(consensustables[["DCIS-normal"]]$curated[ consensustables[["DCIS-normal"]]$consensus])

print( xtable( rbind(nonexpressed,expressed,significant), caption="", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> purple </th> <th> rejected </th> <th> white </th> <th> yellow </th>  </tr>
  <tr> <td align="right"> nonexpressed </td> <td align="right">   70 </td> <td align="right">  123 </td> <td align="right">  356 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> expressed </td> <td align="right">    3 </td> <td align="right">   24 </td> <td align="right">  234 </td> <td align="right">    1 </td> </tr>
  <tr> <td align="right"> significant </td> <td align="right">    1 </td> <td align="right">    9 </td> <td align="right">   68 </td> <td align="right">    1 </td> </tr>
   </table>

In the above table, *nonexpressed* is the number of micorRNAs that were filtered out in one or both of the platforms. The filter was mainly based on some expression threshold in a minimum amount of samples. For short these can be considered non-expressed. The *expressed* class consists of the microRNAs that passed the filter criteria in both platforms. The *significant* class is here the microRNAs that were found as significant in the "normal to DCIS"-comparison for the combined METABRIC and AHUS data sets.

<br/>
<br/>
<br/>

## Stratification based on iClust-subtypes

The UCAM samples are classified as an iClust-sub-type based on mRNA samples. This is provided as a sample annotation like the pam50 estimates. We will find differentially expressed microRNA between the DCIS samples and each of the iClust-types. This classification is only available for the UCAM data set and no alternative meta-analysis is done.




```r
iClust = sampleannotation$iClust
iClust[!is.na(iClust)] = paste("iClust", iClust[!is.na(iClust)], sep="")
iClust[sampleannotation$tissue_type=="DCIS"] = "DCIS"
to_be_used =  sampleannotation$tissue_type %in% c("DCIS", "invasive") &
	!is.na(iClust) & # do not use samples without clinical info
	sampleannotation$provider=="UCAM" 
table(iClust[to_be_used],sampleannotation[to_be_used,"tissue_type"], useNA="ifany")
```

```
##           
##            DCIS invasive
##   DCIS       10        0
##   iClust1     0       91
##   iClust10    0      134
##   iClust2     0       52
##   iClust3     0      191
##   iClust4     0      242
##   iClust5     0      122
##   iClust6     0       56
##   iClust7     0      114
##   iClust8     0      186
##   iClust9     0       87
```


The diff test, result files saved to [here](output-analysis/used_in_article/additionalfile1/difflists_iClust_UCAM_only).

```r
labels = unique(iClust[to_be_used & !iClust %in% c("", "unknown", "DCIS")])
contrasts = paste(labels, "-DCIS", sep="")

group = factor(iClust[to_be_used])

design = model.matrix(~0+group)
colnames(design) = gsub("group", "", colnames(design))
fit = lmFit(common_matrix[,to_be_used], design)
fit$genes$name = miRNA2MIMAT[rownames(common_matrix), "preferredname"]
fit$genes$MIMAT = rownames(common_matrix)
fit$genes$curated = curatedtab[ match( rownames(common_matrix), curatedtab$MIMAT ), "status"]
cont.matrix = makeContrasts ( contrasts=contrasts, levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

iclustres = list()
for(thiscontrast in contrasts)
{
  thisres = topTable(fit2, coef=thiscontrast, adjust="BH", number=999)

  names(thisres)=gsub("^ID.", "", names(thisres))
  thisname= gsub("-", "-vs-", thiscontrast)
  thislabel = strsplit(thiscontrast, "-")[[1]][1]
	
  iclustres[[thislabel]] = thisres
  
  difflistfile=paste(iClustdifflistdir, "/difflist-iClust-",thisname, postfix, sep="")
  write.table(thisres, 
            file=difflistfile,
            quote=FALSE, sep=separator, row.names=FALSE, col.names=TRUE)
  print( paste("Diff genes with fdr<0.05 for ", thislabel, " vs DCIS: ",
               sum(thisres$adj.P.Val < 0.05), sep=""))
  
}
```

```
## [1] "Diff genes with fdr<0.05 for iClust10 vs DCIS: 73"
## [1] "Diff genes with fdr<0.05 for iClust4 vs DCIS: 0"
## [1] "Diff genes with fdr<0.05 for iClust2 vs DCIS: 9"
## [1] "Diff genes with fdr<0.05 for iClust1 vs DCIS: 20"
## [1] "Diff genes with fdr<0.05 for iClust3 vs DCIS: 0"
## [1] "Diff genes with fdr<0.05 for iClust5 vs DCIS: 1"
## [1] "Diff genes with fdr<0.05 for iClust9 vs DCIS: 11"
## [1] "Diff genes with fdr<0.05 for iClust6 vs DCIS: 10"
## [1] "Diff genes with fdr<0.05 for iClust7 vs DCIS: 1"
## [1] "Diff genes with fdr<0.05 for iClust8 vs DCIS: 17"
```
## Printing more tables and plots used in the article.

Creating table 7 from our article. All the microRNAs that are found to change both from normal to DCIS and further from DCIS to one of the sub-types of invasive. Also creating Figure 3 which highlights 4 microRNA.


```r
dcisDEG = as.character(subset(volinia_overlap[["DCIS-normal"]]$MIMATID_vol, volinia_overlap[["DCIS-normal"]]$validated))

table7 = data.frame(MIMA=character(0), name=character(0),  DCIS_direction=character(0),
										subtype=character(0), annottype=character(0), subtype_direction=character(0),
										stringsAsFactors=FALSE)


for(annot in c("pam50_est", "IHC"))
{
	for(lab in unique(sampleannotation[,annot]))
	{

		thislab = paste(lab, "-DCIS", sep="")
		a = consensustables[[thislab]]$consensus == TRUE &
				consensustables[[thislab]]$curated %in% c("purple", "white", "yellow")
		a = a & !is.na(a)
		commonDEG = intersect(consensustables[[thislab]]$MIMAT[a], dcisDEG)
		if(length(commonDEG) > 0)
		{
			df = data.frame( MIMAT = commonDEG,
										name=miRNA2MIMAT[commonDEG, "preferredname"],
										DCIS_direction=ifelse(consensustables[["DCIS-normal"]][commonDEG, "logFC.merged"] > 0, "UP", "DOWN" ),
										annottype=annot,
										subtype=lab,
										subtype_direction=ifelse(consensustables[[thislab]][commonDEG, "logFC.merged"] > 0, "UP", "DOWN" ),
										stringsAsFactors=FALSE
										)		
			table7=rbind(table7,df )
		}

	}
}

for(lab in names(iclustres))
{

		a = iclustres[[lab]]$adj.P.Val < fdrCO & iclustres[[lab]]$curated %in% c("purple", "white", "yellow")
		commonDEG = intersect( iclustres[[lab]]$MIMAT[a], dcisDEG)
		a = a & !is.na(a)
		if(length(commonDEG) > 0)
		{
			df = data.frame( MIMAT = commonDEG,
										name=miRNA2MIMAT[commonDEG, "preferredname"],
										DCIS_direction=ifelse(consensustables[["DCIS-normal"]][commonDEG, "logFC.merged"] > 0, "UP", "DOWN" ),
										annottype="iClust",
										subtype=lab,
										subtype_direction=ifelse(iclustres[[lab]][commonDEG, "logFC"] > 0, "UP", "DOWN" ),
										stringsAsFactors=FALSE
										)		
			table7=rbind(table7,df )
		}
}

nicetable7 = aggregate(subtype ~., table7[, -4], FUN=paste, collapse=" ")
nicetable7 = nicetable7[order(nicetable7$name),]
nicetable7 = nicetable7[order(nicetable7$DCIS_direction, decreasing=TRUE),]

write.table(nicetable7, paste(articledir, "/table7", postfix, sep=""), quote=TRUE, sep=separator, row.names=FALSE)

pdf(file=paste(supplementaryfile1dir, "/boxplots_table7microRNA.pdf", sep=""))

par(oma=c(0,0,2,0))
par(mfrow=c(2, 1))
par(mar=c(2, 2, 2, 0.5) + 0.1)
par(mgp = c(1.2, 0.2, 0))
boxplotcex = 1.25
for(thismimat in unique(table7$MIMAT))
{
	tmptab = 	table7[table7$MIMAT==thismimat,]
	
	for(p in unique(sampleannotation$provider))
	{
		datapoints = list()		
		datapoints[["normal"]] = common_matrix[thismimat, sampleannotation$tissue_type=="normal" & sampleannotation$provider==p ]
		datapoints[["DCIS"]] = common_matrix[thismimat, sampleannotation$tissue_type=="DCIS" & sampleannotation$provider==p ]
		boxcols = c("white", "grey85")
		
		for(i in 1:nrow(tmptab))
		{
			n = sampleannotation[,tmptab$annottype[i] ] == gsub("iClust", "", tmptab$subtype[i]) & sampleannotation$provider==p
			n = n & !is.na(n)
			if(sum(n)>0)
			{
				datapoints[[tmptab$subtype[i]]] = common_matrix[thismimat, n]
				boxcols = c(boxcols, "grey45")
			}																					 
		}
		
		boxplot(datapoints, ylab="Log2 processed signal", col=boxcols,  cex.axis=boxplotcex, cex.lab=0.7, yaxt="n", xaxt="n")
		axis(1, at=1:length(datapoints), labels=names(datapoints), cex.axis=boxplotcex, tick=FALSE)
		axis(2, cex.axis=0.7, tick=FALSE)
		legend("topright", legend=p, cex=boxplotcex, bty="n")
	}
	title(main=paste(thismimat, miRNA2MIMAT[thismimat, "preferredname"], sep=",   "),outer=T)
}
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Figure3, a subset of 4 microRNAS
pdf(file=paste(articledir, "/fig3.pdf", sep=""))

par(oma=c(0,0,0,0))
par(mfcol=c(4, 2))
par(mgp = c(1.2, 0.2, 0))
boxplotcex = 1
for(thismicroRNA in c("hsa-miR-193a-5p", "hsa-miR-378a-3p", "hsa-miR-342-3p", "hsa-miR-497-5p"))
{
	tmptab = 	table7[table7$name==thismicroRNA,]
	thismimat = tmptab$MIMAT[1]
	nicename = gsub("hsa-miR", "mir", thismicroRNA)
	
	for(p in unique(sampleannotation$provider))
	{
		datapoints = list()		
		datapoints[["normal"]] = common_matrix[thismimat, sampleannotation$tissue_type=="normal" & sampleannotation$provider==p ]
		datapoints[["DCIS"]] = common_matrix[thismimat, sampleannotation$tissue_type=="DCIS" & sampleannotation$provider==p ]
		boxcols = c("white", "grey85")
		
		for(i in 1:nrow(tmptab))
		{
			n = sampleannotation[,tmptab$annottype[i] ] == gsub("iClust", "", tmptab$subtype[i]) & sampleannotation$provider==p
			n = n & !is.na(n)
			if(sum(n)>0)
			{
				datapoints[[tmptab$subtype[i]]] = common_matrix[thismimat, n]
				boxcols = c(boxcols, "grey45")
			}																					 
		}
		if(p=="AHUS")
		{
				par(mar=c(1, 1, 2, 0.5) + 0.1)
				plottitle=nicename
		}else{
				par(mar=c(2, 1, 1, 0.5) + 0.1)
				plottitle=""
		}
		boxplot(datapoints, col=boxcols,  cex.axis=boxplotcex, cex.lab=0.7, yaxt="n", xaxt="n", main=plottitle, cex.main=2)
		axis(1, at=1:length(datapoints), labels=names(datapoints), cex.axis=boxplotcex, tick=FALSE)
		axis(2, cex.axis=0.7, tick=FALSE)
		legend("topright", legend=p, cex=1.5, bty="n")
	}
	#title(main=paste(thismimat, miRNA2MIMAT[thismimat, "preferredname"], sep=",   "),outer=T)
}
dev.off()
```

```
## quartz_off_screen 
##                 2
```



```r
sessionInfo()
```

```
R version 3.1.1 (2014-07-10)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] MAMA_2.2.1            GeneMeta_1.38.0       gtools_3.4.2         
 [4] multtest_2.22.0       metaMA_2.1            SMVar_1.3.3          
 [7] sva_3.12.0            genefilter_1.48.1     mgcv_1.8-6           
[10] nlme_3.1-120          xtable_1.7-4          RColorBrewer_1.1-2   
[13] knitr_1.10.5          AgiMicroRna_2.16.0    affycoretools_1.38.0 
[16] GO.db_3.0.0           RSQLite_1.0.0         DBI_0.3.1            
[19] AnnotationDbi_1.28.2  GenomeInfoDb_1.2.5    IRanges_2.0.1        
[22] S4Vectors_0.4.0       preprocessCore_1.28.0 affy_1.44.0          
[25] limma_3.22.7          Biobase_2.26.0        BiocGenerics_0.12.1  
[28] plyr_1.8.2           

loaded via a namespace (and not attached):
 [1] acepack_1.3-3.3           affyio_1.34.0            
 [3] annaffy_1.38.0            annotate_1.44.0          
 [5] AnnotationForge_1.8.2     base64enc_0.1-2          
 [7] BatchJobs_1.6             BBmisc_1.9               
 [9] BiocInstaller_1.16.5      BiocParallel_1.0.3       
[11] biomaRt_2.22.0            Biostrings_2.34.1        
[13] biovizBase_1.14.1         bit_1.1-12               
[15] bitops_1.0-6              brew_1.0-6               
[17] BSgenome_1.34.1           Category_2.32.0          
[19] caTools_1.17.1            checkmate_1.5.3          
[21] cluster_2.0.1             codetools_0.2-11         
[23] colorspace_1.2-6          DESeq2_1.6.3             
[25] dichromat_2.0-0           digest_0.6.8             
[27] edgeR_3.8.6               evaluate_0.7             
[29] fail_1.2                  ff_2.2-13                
[31] foreach_1.4.2             foreign_0.8-63           
[33] formatR_1.2               Formula_1.2-1            
[35] gcrma_2.38.0              gdata_2.16.1             
[37] geneplotter_1.44.0        GenomicAlignments_1.2.2  
[39] GenomicFeatures_1.18.7    GenomicRanges_1.18.4     
[41] GGally_0.5.0              ggbio_1.14.0             
[43] ggplot2_1.0.1             GOstats_2.32.0           
[45] gplots_2.17.0             graph_1.44.1             
[47] gridExtra_0.9.1           GSEABase_1.28.0          
[49] gtable_0.1.2              Hmisc_3.16-0             
[51] hwriter_1.3.2             iterators_1.0.7          
[53] KernSmooth_2.23-14        lattice_0.20-31          
[55] latticeExtra_0.6-26       locfit_1.5-9.1           
[57] magrittr_1.5              markdown_0.7.7           
[59] MASS_7.3-40               Matrix_1.2-0             
[61] MergeMaid_2.38.0          metaArray_1.44.0         
[63] mime_0.3                  munsell_0.4.2            
[65] nnet_7.3-9                oligoClasses_1.28.0      
[67] OrganismDbi_1.8.1         PFAM.db_3.0.0            
[69] proto_0.3-10              R.methodsS3_1.7.0        
[71] R.oo_1.19.0               R.utils_2.0.2            
[73] RBGL_1.42.0               Rcpp_0.11.6              
[75] RcppArmadillo_0.5.100.1.0 RCurl_1.95-4.6           
[77] ReportingTools_2.6.0      reshape_0.8.5            
[79] reshape2_1.4.1            rpart_4.1-9              
[81] Rsamtools_1.18.3          rtracklayer_1.26.3       
[83] scales_0.2.4              sendmailR_1.2-1          
[85] splines_3.1.1             stringi_0.4-1            
[87] stringr_1.0.0             survival_2.38-1          
[89] tools_3.1.1               VariantAnnotation_1.12.9 
[91] XML_3.98-1.1              XVector_0.6.0            
[93] zlibbioc_1.12.0          
```

generation ended 2015-10-28 12:45:27. Time spent 0 minutes .



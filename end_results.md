Overlap between our results and Volinia et al.
========================================================
2015-06-10 14:08:54


<br/>
<br/>

We will check for

- concordance between our two approaches called **meta** and **merged**
- concordance between our results and Volinia et al.
- microRNA evidence status for the results

<BR/>
<BR/>
Read in prepared data.


```r
library(xtable)
if(!exists("inputisread"))
	source("read_input.r")
postfix = ".txt"
options(digits=4)
outputdir = "output-analysis"
consensusdir = paste(outputdir, "/consensus", sep="")
if(!file.exists(consensusdir))
	dir.create(consensusdir)
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
<!-- Wed Jun 10 14:08:54 2015 -->
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
<!-- Wed Jun 10 14:08:54 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> Basallike </th> <th> DCIS </th> <th> Her2 </th> <th> LumA </th> <th> LumB </th> <th> Normallike </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    5 </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">   16 </td> <td align="right">   15 </td> <td align="right">    7 </td> <td align="right">    4 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun 10 14:08:54 2015 -->
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
<!-- Wed Jun 10 14:08:54 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> HER2neg_ERneg_PGRneg </th> <th> HER2neg_ERpos </th> <th> HER2pos_ERneg </th> <th> HER2pos_ERpos </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> benign </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   23 </td> </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    8 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    1 </td> <td align="right">   28 </td> <td align="right">    4 </td> <td align="right">    8 </td> <td align="right">   14 </td> </tr>
  <tr> <td align="right"> normal </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   70 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun 10 14:08:54 2015 -->
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
write.table(table1, file=paste(outputdir, "/table1.txt", sep=""), sep="\t", col.names=NA, quote=FALSE)
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
  thistable = cbind(MIMAT = thistable$Row.names, name=unfilteredcommonname[thistable$Row.names], thistable[,-1])
  thistable$same_dir =  (thistable$logFC.merged * thistable$TestStatistic.meta ) > 0
  thistable$consensus = thistable$adj.P.Val.merged <= 0.05 & 
  	thistable$adj.P.Val.meta <= 0.05 & thistable$same_dir
  consensustables[[thiscomp]] = thistable[order(thistable$consensus, -thistable$adj.P.Val.merged, decreasing=TRUE),]
}

if(sum(summarytab$wrong_direction)>0)
	warning("Some common microRNAs were found that did not have the same dirrection of change!!")
```
The above tables are written to [files](output-analysis/consensus) later.
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
<!-- Wed Jun 10 14:08:54 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> Overlap between meta and merged </caption>
<tr> <th>  </th> <th> found_meta </th> <th> found_merged </th> <th> found_both </th>  </tr>
  <tr> <td align="right"> DCIS-normal </td> <td align="right">  116 </td> <td align="right">   86 </td> <td align="right">   82 </td> </tr>
  <tr> <td align="right"> invasive-DCIS </td> <td align="right">    0 </td> <td align="right">   10 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> Normallike-DCIS </td> <td align="right">    7 </td> <td align="right">   34 </td> <td align="right">    7 </td> </tr>
  <tr> <td align="right"> LumA-DCIS </td> <td align="right">   10 </td> <td align="right">   31 </td> <td align="right">    7 </td> </tr>
  <tr> <td align="right"> LumB-DCIS </td> <td align="right">   30 </td> <td align="right">   43 </td> <td align="right">   24 </td> </tr>
  <tr> <td align="right"> Her2-DCIS </td> <td align="right">   11 </td> <td align="right">   16 </td> <td align="right">    6 </td> </tr>
  <tr> <td align="right"> Basallike-DCIS </td> <td align="right">    7 </td> <td align="right">   61 </td> <td align="right">    7 </td> </tr>
  <tr> <td align="right"> HER2neg_ERneg_PGRneg-DCIS </td> <td align="right">    2 </td> <td align="right">   30 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> HER2neg_ERpos-DCIS </td> <td align="right">    1 </td> <td align="right">   11 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> HER2pos_ERneg-DCIS </td> <td align="right">   15 </td> <td align="right">   19 </td> <td align="right">    8 </td> </tr>
  <tr> <td align="right"> HER2pos_ERpos-DCIS </td> <td align="right">    3 </td> <td align="right">    5 </td> <td align="right">    0 </td> </tr>
   </table>


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
signchangetab[,"name"] = unfilteredcommonname[consensustables[[1]]$MIMAT]

for(n in names(consensustables))
{
	x = consensustables[[n]][rownames(signchangetab),]
	signchangetab[ x$consensus & x$logFC.merged<0, n ] = "DOWN"
	signchangetab[ x$consensus & x$logFC.merged>0, n ] = "UP" 
}
signchangetabfn = paste(consensusdir, "/signchangetab", postfix, sep="")
write.table(signchangetab, file=signchangetabfn,
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="")
```

[Link to file](output-analysis/consensus/signchangetab.txt)
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

correlationmatrixfn = paste(consensusdir, "/subtypecorrelation", postfix, sep="")
h1 = "Count of significant microRNAs found going in the same direction for the subtypes"
write(h1, file=correlationmatrixfn)
write(paste("\t", colnames(correlationmatrix1), collapse="\t"), file=correlationmatrixfn)
write.table(correlationmatrix1, file=correlationmatrixfn,quote=FALSE, sep="\t", col.names=FALSE, na="X", append=TRUE)
write("\n\n\n", file=correlationmatrixfn, append=TRUE)
h2 = "Count of significant microRNAs found going in the opposite direction for the subtypes"
write(h2, file=correlationmatrixfn, append=TRUE)
write(paste("\t", colnames(correlationmatrix2), collapse="\t"), file=correlationmatrixfn)
write.table(correlationmatrix2, file=correlationmatrixfn,quote=FALSE, sep="\t", col.names=FALSE,  na="X", append=TRUE)

print(xtable(correlationmatrix1, 
             caption=h1, digits=3), 
      comment = TRUE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun 10 14:08:54 2015 -->
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

[Link to file](output-analysis/consensus/subtypecorrelation.txt)

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
                  list(MIMATID=mimatmapping[paste("hsa-", voliniatab[[n]]$miRNA, sep="")]), 
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
##  miR-218(2)         <NA> -0.6897 1.59e-05
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
            quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)

	rm(x, overlaptab, contrast)
}

# write table2 in our article
table2 = consensustables[["DCIS-normal"]]
table2 = table2[table2$consensus,]
table2 = table2[table2$MIMAT %in% rownames( volinia_overlap[["DCIS-normal"]])[volinia_overlap[["DCIS-normal"]]$validated], ]
table2 = table2[order(table2$logFC.merged),2:6]
table2[,c(2,4)] = round(table2[,c(2,4)], 2)
write.table(table2, file=paste(outputdir, "/table2.txt", sep=""), sep="\t", col.names=NA, quote=FALSE)
```

And below is the result list from Volinia et al. with our results attached. First the "DCIS-normal" comparison. Rows are sorted by FDR from our "merged" approach. [Link to file](output-analysis/Volinia_validation_DCIS-normal.txt)


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
  <tr> <td> miR-378 </td> <td> MIMAT0000731_MIMAT0000732 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-155 </td> <td> MIMAT0000646 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-200c </td> <td> MIMAT0000617 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-425 </td> <td> MIMAT0001343_MIMAT0003393 </td> <td> UP </td> <td align="right"> 0.012 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
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
  <tr> <td> miR-106b </td> <td> MIMAT0000680 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-145* </td> <td> MIMAT0004601 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-342 </td> <td> MIMAT0000753 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-100 </td> <td> MIMAT0000098 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.008 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> DOWN </td> <td align="right"> 0.009 </td> <td> UP </td> <td align="right"> 0.022 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-22 </td> <td> MIMAT0000077 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.024 </td> <td> DOWN </td> <td align="right"> 0.049 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-652 </td> <td> MIMAT0003322 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-374b </td> <td> MIMAT0004955 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.051 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-19a </td> <td> MIMAT0000073 </td> <td> UP </td> <td align="right"> 0.008 </td> <td> UP </td> <td align="right"> 0.073 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-376c </td> <td> MIMAT0000720 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> DOWN </td> <td align="right"> 0.088 </td> <td> DOWN </td> <td align="right"> 0.324 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.112 </td> <td> UP </td> <td align="right"> 0.092 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-361-5p </td> <td> MIMAT0000703 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.114 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-30d </td> <td> MIMAT0000245 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.116 </td> <td> DOWN </td> <td align="right"> 0.034 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-223 </td> <td> MIMAT0000280 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> UP </td> <td align="right"> 0.124 </td> <td> UP </td> <td align="right"> 0.079 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-185 </td> <td> MIMAT0000455 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.130 </td> <td> UP </td> <td align="right"> 0.183 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-374a </td> <td> MIMAT0000727 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.195 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-127-3p </td> <td> MIMAT0000446 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.200 </td> <td> DOWN </td> <td align="right"> 0.084 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-26b </td> <td> MIMAT0000083 </td> <td> UP </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.377 </td> <td> UP </td> <td align="right"> 0.265 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-574-3p </td> <td> MIMAT0003239 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.378 </td> <td> UP </td> <td align="right"> 0.271 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-20a </td> <td> MIMAT0000075 </td> <td> UP </td> <td align="right"> 0.007 </td> <td> UP </td> <td align="right"> 0.446 </td> <td> UP </td> <td align="right"> 0.350 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-107 </td> <td> MIMAT0000104 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.468 </td> <td> UP </td> <td align="right"> 0.583 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-125a </td> <td> MIMAT0000443 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.782 </td> <td> UP </td> <td align="right"> 0.628 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-423-5p </td> <td> MIMAT0004748 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.852 </td> <td> UP </td> <td align="right"> 0.756 </td> <td> FALSE </td> </tr>
  <tr> <td> let-7d </td> <td> MIMAT0000065 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.887 </td> <td> DOWN </td> <td align="right"> 0.731 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15b </td> <td> MIMAT0000417 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.945 </td> <td> UP </td> <td align="right"> 0.641 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-29c </td> <td> MIMAT0000681 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.951 </td> <td> DOWN </td> <td align="right"> 0.639 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15a </td> <td> MIMAT0000068 </td> <td> UP </td> <td align="right"> 0.013 </td> <td> DOWN </td> <td align="right"> 0.972 </td> <td> DOWN </td> <td align="right"> 0.981 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-182 </td> <td> MIMAT0000259 </td> <td> UP </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-423-3p </td> <td> MIMAT0001340 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-106b* </td> <td> MIMAT0004672 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-32 </td> <td> MIMAT0000090 </td> <td> UP </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-28-3p </td> <td> MIMAT0004502 </td> <td> DOWN </td> <td align="right"> 0.004 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-340 </td> <td> MIMAT0000750_MIMAT0004692 </td> <td> UP </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-452 </td> <td> MIMAT0001635 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-17* </td> <td> MIMAT0000071 </td> <td> DOWN </td> <td align="right"> 0.015 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
   </table>

<br/>
<br/>
And here is the same table for the "invasive-DCIS" comparison. [Link to file](output-analysis/Volinia_validation_invasive-DCIS.txt)


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
  <tr> <td> let-7d </td> <td> MIMAT0000065 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.057 </td> <td> UP </td> <td align="right"> 0.713 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.300 </td> <td> UP </td> <td align="right"> 0.932 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.307 </td> <td> DOWN </td> <td align="right"> 0.713 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-126 </td> <td> MIMAT0000445 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.361 </td> <td> DOWN </td> <td align="right"> 0.824 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-143 </td> <td> MIMAT0000435 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.499 </td> <td> DOWN </td> <td align="right"> 0.958 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-335-5p </td> <td> MIMAT0000765 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.616 </td> <td> DOWN </td> <td align="right"> 0.958 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-10b </td> <td> MIMAT0000254 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.991 </td> <td> UP </td> <td align="right"> 0.958 </td> <td> FALSE </td> </tr>
   </table>
<br/>
<br/>


<br/>
<br/>

## Checking results towards a curated list of microRNA´s

Bastian Fromm and colleagues have carefully evaluated the evidence for the Human micro RNAs reported in mirBase and provided a classification based on this. In short, a lot of microRNAs in mirBASE and included in the microarray platforms are probably not real microRNAs.

A list of human microRNAs and their curation classification is provided:


```r
	curatedfn = "curated_microRNA_MIMAT.csv"
  orgcuratedtab  = read.table(paste(annotdir, "/", curatedfn, sep=""), sep="\t", 
  												 header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="")
	orgcuratedtab$miRNA = paste("hsa-", orgcuratedtab$miRNA, sep="")
  # ad hoc fix
	#orgcuratedtab$miRNA = gsub("-mir-", "-miR-", orgcuratedtab$miRNA)
	
	curatedtab = orgcuratedtab[,c(1,2,3)]
	names(curatedtab)[3]= "MIMAT"
	
	curatedtab = rbind(curatedtab, 
										 setNames( orgcuratedtab[,c(1,2,4)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,5)], names(curatedtab)),
										 setNames( orgcuratedtab[,c(1,2,6)], names(curatedtab)) )
	curatedtab = curatedtab[curatedtab$MIMAT != "", ]

	#curatedtab$MIMAT = mimatmapping[curatedtab$miRNA]
	#table(curatedtab$status)
	#curatedtab = curatedtab[!is.na(curatedtab$MIMAT),]
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
## [1] "invasive-DCIS white 7"
## [1] "DCIS-normal white 55"
```
The microRNAs reported by Volinia et al. are all in the "for sure microRNA" list in Fromm et al.
<br/>
<br/>
And now the same assignment and summary for our results from the AHUS and METABRIC data sets.


```r
namelists = list()
for(n in names(consensustables))
{
	x = unlist(lapply(rownames(consensustables[[n]]), FUN=function(x)strsplit(x, "_")[[1]][1]))
	consensustables[[n]]$curated = curatedtab[ match( x, curatedtab$MIMAT ), "status"]
	write.table( format( consensustables[[n]], digits=4), 
            file=paste(consensusdir, "/consensusresults_",meta2merged[n, "comparisontype"], "__", n, postfix, sep=""),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	
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

# write table3 and table4 in our article
tableclasses = c("LumA-DCIS", "LumB-DCIS", "Her2-DCIS", "Basallike-DCIS", "Normallike-DCIS")
nrows = max( unlist(lapply( namelists[tableclasses], FUN=length) ))
table3 = matrix("",ncol=length(tableclasses), nrow=nrows)
colnames(table3)=tableclasses
for(n in tableclasses)
	table3[1:length(namelists[[n]]), n] = namelists[[n]]
write.table(table3, file=paste(outputdir, "/table3.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

tableclasses = c("HER2pos_ERpos-DCIS", "HER2neg_ERpos-DCIS", "HER2pos_ERneg-DCIS", "HER2neg_ERneg_PGRneg-DCIS")
nrows = max( unlist(lapply( namelists[tableclasses], FUN=length) ))
table4 = matrix("",ncol=length(tableclasses), nrow=nrows)
colnames(table4)=tableclasses
for(n in tableclasses)
	if(length(namelists[[n]])>0)
		table4[1:length(namelists[[n]]), n] = namelists[[n]]
write.table(table4, file=paste(outputdir, "/table4.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
```

Result files for all the comparisons with p-values for the two approaches, and with the curation status added, are [here](output-analysis/consensus). The tables look like this (top of the DCIS-normal comparison as an example):


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
  <tr> <td align="right"> MIMAT0005893 </td> <td> MIMAT0005893 </td> <td> hsa-miR-1305 </td> <td align="right"> 0.681 </td> <td align="right"> 0.000 </td> <td align="right"> 6.121 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> purple </td> </tr>
  <tr> <td align="right"> MIMAT0005942 </td> <td> MIMAT0005942 </td> <td> hsa-miR-1288 </td> <td align="right"> 0.342 </td> <td align="right"> 0.000 </td> <td align="right"> 5.097 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> rejected </td> </tr>
  <tr> <td align="right"> MIMAT0000076 </td> <td> MIMAT0000076 </td> <td> hsa-miR-21 </td> <td align="right"> 2.407 </td> <td align="right"> 0.000 </td> <td align="right"> 8.676 </td> <td align="right"> 0.000 </td> <td> TRUE </td> <td> TRUE </td> <td> white </td> </tr>
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
  <tr> <td align="right"> Normallike-DCIS </td> <td align="right">    4 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> LumA-DCIS </td> <td align="right">    6 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    1 </td> <td align="right">    0 </td> </tr>
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
filt1 = unlist(lapply(filteredcommonMIMAT, FUN=function(x)strsplit(x, "_")[[1]][1]))
unfilt1 = unlist(lapply(unfilteredcommonMIMAT, FUN=function(x)strsplit(x, "_")[[1]][1]))

expressed =	table(curatedtab$status[match(filt1, curatedtab$MIMAT)])
nonexpressed=	table(curatedtab$status[match(unfilt1[!unfilt1 %in% filt1],curatedtab$MIMAT) ])

#genelist = rownames(consensustables[["DCIS-normal"]])[consensustables[["DCIS-normal"]]$consensus]
#significant = table(curatedtab$status[curatedtab$MIMAT %in% genelist])

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
  <tr> <td align="right"> nonexpressed </td> <td align="right">   70 </td> <td align="right">  122 </td> <td align="right">  350 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> expressed </td> <td align="right">    3 </td> <td align="right">   24 </td> <td align="right">  235 </td> <td align="right">    1 </td> </tr>
  <tr> <td align="right"> significant </td> <td align="right">    1 </td> <td align="right">    9 </td> <td align="right">   68 </td> <td align="right">    1 </td> </tr>
   </table>

In the above table, *nonexpressed* is the number of micorRNAs that were filtered out in one or both of the platforms. The filter was mainly based on some expression threshold in a minimum amount of samples. For short these can be considered non-expressed. The *expressed* class consists of the microRNAs that passed the filter criteria in both platforms. The *significant* class is here the microRNAs that were found as significant in the "normal to DCIS"-comparison for the combined METABRIC and AHUS data sets.

<br/>
<br/>


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
[13] RColorBrewer_1.1-2    knitr_1.9             AgiMicroRna_2.14.0   
[16] affycoretools_1.36.1  GO.db_2.14.0          RSQLite_1.0.0        
[19] DBI_0.3.1             AnnotationDbi_1.26.1  GenomeInfoDb_1.0.2   
[22] preprocessCore_1.26.1 affy_1.42.3           limma_3.20.9         
[25] Biobase_2.24.0        BiocGenerics_0.10.0   data.table_1.9.2     
[28] plyr_1.8.1           

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
[69] R.oo_1.19.0               R.utils_2.0.2            
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

generation ended 2015-06-10 14:08:55. Time spent 0 minutes .



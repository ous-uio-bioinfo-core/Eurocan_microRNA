Overlap between our results and Volinia et al.
========================================================
2015-06-03 17:45:02


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
<!-- Wed Jun  3 17:45:03 2015 -->
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
<!-- Wed Jun  3 17:45:03 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> Basallike </th> <th> DCIS </th> <th> Her2 </th> <th> LumA </th> <th> LumB </th> <th> Normallike </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    5 </td> <td align="right">    0 </td> <td align="right">    8 </td> <td align="right">   16 </td> <td align="right">   15 </td> <td align="right">    7 </td> <td align="right">    4 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun  3 17:45:03 2015 -->
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
<!-- Wed Jun  3 17:45:03 2015 -->
<table CELLPADDING=5>
<caption align="bottom"> AHUS </caption>
<tr> <th>  </th> <th> HER2neg_ERneg_PGRneg </th> <th> HER2neg_ERpos </th> <th> HER2pos_ERneg </th> <th> HER2pos_ERpos </th> <th> unknown </th>  </tr>
  <tr> <td align="right"> benign </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   23 </td> </tr>
  <tr> <td align="right"> DCIS </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    8 </td> </tr>
  <tr> <td align="right"> invasive </td> <td align="right">    1 </td> <td align="right">   28 </td> <td align="right">    4 </td> <td align="right">    8 </td> <td align="right">   14 </td> </tr>
  <tr> <td align="right"> normal </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">    0 </td> <td align="right">   70 </td> </tr>
   </table>
<!-- html table generated in R 3.1.1 by xtable 1.7-4 package -->
<!-- Wed Jun  3 17:45:03 2015 -->
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




































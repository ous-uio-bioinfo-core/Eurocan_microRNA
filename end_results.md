Overlap between data sets and Volinia.
========================================================
2014-12-10 16:24:47


<br/>
<br/>

Agreement between our two approaches called **merged** and **meta**.
Agreement between our results and Volinia.
Checking all results toward a list of verified microRNAs.

<BR/>
<BR/>
## Meta vs Merged results.

Our data sets came from two providers, "AHUS" and "UCAM", using two different Agilent microRNA platforms. Batch effects were present as visualized in the QC report. We tested for differentially expressed genes using two approaches, combining all data in Limma using provider as a blocking factor ("merged" approach), and analyzing the two sets inependetly with Limma and combine the p-values ("meta" approach).

Here we will list the agreemet between the two apporaches.



```r
#metafilefolder = "not_in_github/t-test"
metafilefolder = "not_in_github/lilianas"
mergedfilefolder = "output-analysis"
#metaresultfiles = list.files("lilianas/1000_permutations")
#mergedresultfiles = list.files("output-analysis-2014-08-05", pattern="difflist")

comparisons = c("DCIS-normal", "invasive-DCIS", "Normallike-DCIS", "LumA-DCIS", "LumB-DCIS", 
								"Her2-DCIS", "Basallike-DCIS")

# meta2merged = data.frame(metafn=character(0), mergedfn=character(0), direction=numeric(0), stringsAsFactors=FALSE)
# meta2merged["DCIS-normal",] = list("DCIS_to_normal_test", "difflist-DCIS-vs-normal.txt",-1)
# meta2merged["invasive-DCIS",] = list("DCIS_to_invasive_test", "difflist-invasive-vs-DCIS.txt",1)
# meta2merged["Normallike-DCIS",] = list("DCIS_to_normal(invasive)_test", "difflist-pam50-Normallike-vs-DCIS.txt",1)
# meta2merged["LumA-DCIS",] = list("DCIS_to_lumA_test", "difflist-pam50-LumA-vs-DCIS.txt",1)
# meta2merged["LumB-DCIS",] = list("DCIS_to_lumB_test", "difflist-pam50-LumB-vs-DCIS.txt",1)
# meta2merged["Her2-DCIS",] = list("DCIS_to_Her2_test", "difflist-pam50-Her2-vs-DCIS.txt",1)
# meta2merged["Basallike-DCIS",] = list("basal_to_DCIS_test", "difflist-pam50-Basallike-vs-DCIS.txt",-1)


meta2merged = data.frame(metafn=character(0), mergedfn=character(0), direction=numeric(0), stringsAsFactors=FALSE)
meta2merged["DCIS-normal",] = list("DCIS_to_normal_limma.txt", "difflist-DCIS-vs-normal.txt",-1)
meta2merged["invasive-DCIS",] = list("DCIS_to_invasive_limma.txt", "difflist-invasive-vs-DCIS.txt",1)
meta2merged["Normallike-DCIS",] = list("dcis_to_normal_like(invasive)_limma.txt", "difflist-pam50-Normallike-vs-DCIS.txt",1)
meta2merged["LumA-DCIS",] = list("dcis_to_lumA_limma.txt", "difflist-pam50-LumA-vs-DCIS.txt",1)
meta2merged["LumB-DCIS",] = list("dcis_to_lumB_limma.txt", "difflist-pam50-LumB-vs-DCIS.txt",1)
meta2merged["Her2-DCIS",] = list("dcis_to_her2_limma.txt", "difflist-pam50-Her2-vs-DCIS.txt",1)
meta2merged["Basallike-DCIS",] = list("basallike_to_dcis_limma.txt", "difflist-pam50-Basallike-vs-DCIS.txt",-1)

meta2merged["HER2neg_ERneg_PGRneg-DCIS",] = list("dcis_to_HER2neg_ERneg_PGRneg_limma.txt", "difflist-IHC-HER2neg_ERneg_PGRneg-vs-DCIS.txt",1)
meta2merged["HER2neg_ERpos-DCIS",] = list("dcis_to_HER2neg_ERpos_limma.txt", "difflist-IHC-HER2neg_ERpos-vs-DCIS.txt",1)
meta2merged["HER2pos_ERneg-DCIS",] = list("dcis_to_HER2pos_ERneg_limma.txt", "difflist-IHC-HER2pos_ERneg-vs-DCIS.txt",1)
meta2merged["HER2pos_ERpos-DCIS",] = list("dcis_to_HER2pos_ERpos_limma.txt", "difflist-IHC-HER2pos_ERpos-vs-DCIS.txt",1)


metatables = list()
mergedtables = list()
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
    
    #harmonizing
    colnames(thismetatab)[colnames(thismetatab)=="FC"] = "logFC"
    rownames(thismetatab) = thismetatab$MIMAT
    thismetatab$logFC = thismetatab$logFC * meta2merged[thiscomp, "direction"]
  #summarytab[thiscomp, "metafile"] = rownames(thisrows)[1]
  #summarytab[thiscomp, "mergedfile"] = thisrows[1, "mergedfn"]

  summarytab[thiscomp, "found_meta"] = nrow(thismetatab)
  b = thismergedtab$adj.P.Val <= 0.05
  summarytab[thiscomp, "found_merged"] = sum(b)
  
  commonMIMATID = intersect( rownames(thismetatab),  rownames(thismergedtab)[b])
  summarytab[thiscomp, "found_both"] = sum( thismetatab[commonMIMATID, "logFC"] * thismergedtab[commonMIMATID, "logFC"] > 0) 
  rm(commonMIMATID)
  metatables[[thiscomp]] = thismetatab
  mergedtables[[thiscomp]] = thismergedtab
}
```
<br/>
<br/>
<br/>
The number of microRNAs found differentially expressed between states, according to the two approaches. The direction of foldchange is taken into acount. 

```r
print(summarytab[,1:3])
```

```
##                           found_meta found_merged found_both
## DCIS-normal                      113           86         83
## invasive-DCIS                      0           10          0
## Normallike-DCIS                    7           39          7
## LumA-DCIS                         10           31          7
## LumB-DCIS                         33           37         22
## Her2-DCIS                         11           17          6
## Basallike-DCIS                     7           58          7
## HER2neg_ERneg_PGRneg-DCIS          2           31          2
## HER2neg_ERpos-DCIS                 1           13          0
## HER2pos_ERneg-DCIS                15           21          8
## HER2pos_ERpos-DCIS                 3            9          0
```
For several of the comparisons, the shorter list seems to be mostly included in the longer.
<br/>
<br/>


## Our results vs. Volinia et al.´s results


```r
volinia_DCIS_vs_normal_file = "volinia-DCIS-vs-normal.txt"
volinia_invasive_vs_DCIS_file = "volinia-invasive-vs-DCIS.txt"
```

volinia-DCIS-vs-normal.txt and  volinia-invasive-vs-DCIS.txt are copy-pasted from Volinia et al.s [supplementary tables S1 and S2](http://www.pnas.org/content/suppl/2012/02/02/1200010109.DCSupplemental/pnas.201200010SI.pdf). The tables were reformatted from the pdf file into a more practical tabular text format. This is the results we aim to validate with the AHUS and UCAM data sets.


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


```r
print( xtable(head(voliniatab[[1]]), caption="Start of table S2 copied from Volinia et al.", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> Start of table S2 copied from Volinia et al. </caption>
<tr> <th> miRNA </th> <th> MIMATID </th> <th> log2FC </th> <th> FDR </th>  </tr>
  <tr> <td> miR-10b </td> <td> MIMAT0000254 </td> <td align="right"> -1.286 </td> <td align="right"> 0.000 </td> </tr>
  <tr> <td> miR-143 </td> <td> MIMAT0000435 </td> <td align="right"> -1.358 </td> <td align="right"> 0.000 </td> </tr>
  <tr> <td> let-7d </td> <td> MI0000065_MIMAT0000065 </td> <td align="right"> 1.333 </td> <td align="right"> 0.000 </td> </tr>
  <tr> <td> miR-218(2) </td> <td>  </td> <td align="right"> -0.690 </td> <td align="right"> 0.000 </td> </tr>
  <tr> <td> miR-335-5p </td> <td> MIMAT0000765 </td> <td align="right"> -1.184 </td> <td align="right"> 0.000 </td> </tr>
  <tr> <td> miR-126 </td> <td> MIMAT0000445 </td> <td align="right"> -1.120 </td> <td align="right"> 0.000 </td> </tr>
   </table>


Compared to the original tables from Volinia, the fold change is here log2 transformed and "IDC" is renamed to "invasive". 



```r
volinia_overlap = list()

for(contrast in names(voliniatab))
{
	#contrast = "invasive_vs_DCIS"
	overlaptab = voliniatab[[contrast]]
	overlaptab = overlaptab[!is.na(overlaptab$MIMATID),]
	rownames(overlaptab) = overlaptab$MIMATID
	overlaptab$Direction = ifelse(overlaptab$log2FC > 0, "UP", "DOWN")
	overlaptab = overlaptab[, c("miRNA", "MIMATID", "Direction","FDR")]
	names(overlaptab) = paste(names(overlaptab), "_vol", sep="")
		
	tab = mergedtables[[contrast]]
	a = intersect(  rownames(tab), rownames(overlaptab) )
	overlaptab$Direction_merged = NA
	overlaptab[a, "Direction_merged"] = ifelse(tab[a, "logFC"] > 0,  "UP", "DOWN")
	overlaptab$FDR_merged = NA
	overlaptab[a, "FDR_merged"] = tab[a, "adj.P.Val"]

	tab = metatables[[contrast]]
	a = intersect(  rownames(tab), rownames(overlaptab) )
	overlaptab$Direction_meta = NA
	overlaptab[a, "Direction_meta"] = ifelse(tab[a, "logFC"] > 0,  "UP", "DOWN")
	overlaptab$FDR_meta = NA
	overlaptab[a, "FDR_meta"] = fdrCO
	
	x = rowSums (overlaptab[,grepl("Direction", colnames(overlaptab))] == "UP")
	overlaptab$validated = NA
	overlaptab$validated[x %in% c(0,3)] = TRUE
	overlaptab$validated[x %in% c(1,2)] = FALSE
	overlaptab$validated = overlaptab$validated & (overlaptab$FDR_merged<=fdrCO & overlaptab$FDR_meta<=fdrCO)
	
# 	a = overlaptab$Direction_vol > 0 & overlaptab$Direction_merged > 0 & overlaptab$Direction_meta > 0
# 	b = overlaptab$Direction_vol < 0 & overlaptab$Direction_merged < 0 & overlaptab$Direction_meta < 0
# 	c = overlaptab$FDR_merged < 0.05 & overlaptab$FDR_meta < 0.05
# 	overlaptab$validated = (a|b) & c 		
	
	volinia_overlap[[contrast]] = overlaptab
}
```



```r
x = order(volinia_overlap[["DCIS-normal"]]$FDR_merged)
print( xtable( volinia_overlap[["DCIS-normal"]][x,] , caption="microRNA found in Volinia compared to validation data", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> microRNA found in Volinia compared to validation data </caption>
<tr> <th> miRNA_vol </th> <th> MIMATID_vol </th> <th> Direction_vol </th> <th> FDR_vol </th> <th> Direction_merged </th> <th> FDR_merged </th> <th> Direction_meta </th> <th> FDR_meta </th> <th> validated </th>  </tr>
  <tr> <td> miR-21 </td> <td> MIMAT0000076 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-497 </td> <td> MIMAT0002820 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-145 </td> <td> MIMAT0000437 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-193a-5p </td> <td> MIMAT0004614 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-378 </td> <td> MIMAT0000731_MIMAT0000732 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-155 </td> <td> MIMAT0000646 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-200c </td> <td> MIMAT0000617 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-425 </td> <td> MIMAT0001343_MIMAT0003393 </td> <td> UP </td> <td align="right"> 0.012 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-183 </td> <td> MIMAT0000261 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-200b </td> <td> MIMAT0000318 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-429 </td> <td> MIMAT0001536 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-96 </td> <td> MIMAT0000095 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-142-3p </td> <td> MIMAT0000434 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-143* </td> <td> MIMAT0004599 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-140-3p </td> <td> MIMAT0004597 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-142-5p </td> <td> MIMAT0000433 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-99a </td> <td> MIMAT0000097 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> let-7b </td> <td> MI0000063_MIMAT0000063 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> let-7c </td> <td> MI0000064_MIMAT0000064 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-193b </td> <td> MIMAT0002819 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-106b </td> <td> MIMAT0000680 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-145* </td> <td> MIMAT0004601 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-342 </td> <td> MIMAT0000753 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-100 </td> <td> MIMAT0000098 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.008 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> DOWN </td> <td align="right"> 0.009 </td> <td> UP </td> <td align="right"> 0.022 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-22 </td> <td> MIMAT0000077 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.024 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-652 </td> <td> MIMAT0003322 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> </tr>
  <tr> <td> miR-374b </td> <td> MIMAT0004955 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.051 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-19a </td> <td> MIMAT0000073 </td> <td> UP </td> <td align="right"> 0.008 </td> <td> UP </td> <td align="right"> 0.073 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-376c </td> <td> MIMAT0000720 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> DOWN </td> <td align="right"> 0.088 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.112 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-361-5p </td> <td> MIMAT0000703 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.114 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-30d </td> <td> MIMAT0000245 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.116 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-223 </td> <td> MIMAT0000280 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> UP </td> <td align="right"> 0.124 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-185 </td> <td> MIMAT0000455 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.130 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-374a </td> <td> MIMAT0000727 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.195 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-127-3p </td> <td> MIMAT0000446 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.200 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-26b </td> <td> MIMAT0000083 </td> <td> UP </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.377 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-574-3p </td> <td> MIMAT0003239 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.378 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-20a </td> <td> MIMAT0000075 </td> <td> UP </td> <td align="right"> 0.007 </td> <td> UP </td> <td align="right"> 0.446 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-107 </td> <td> MIMAT0000104 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.468 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-125a </td> <td> MIMAT0000443 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.782 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> </tr>
  <tr> <td> miR-423-5p </td> <td> MIMAT0004748 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.852 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> let-7d </td> <td> MI0000065_MIMAT0000065 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.887 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15b </td> <td> MIMAT0000417 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.945 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-29c </td> <td> MIMAT0000681 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.951 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-15a </td> <td> MIMAT0000068 </td> <td> UP </td> <td align="right"> 0.013 </td> <td> DOWN </td> <td align="right"> 0.972 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> miR-182 </td> <td> MIMAT0000259 </td> <td> UP </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-423-3p </td> <td> MIMAT0001340 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-106b* </td> <td> MIMAT0004672 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-32 </td> <td> MIMAT0000090 </td> <td> UP </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-28-3p </td> <td> MIMAT0004502 </td> <td> DOWN </td> <td align="right"> 0.004 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-340 </td> <td> MIMAT0000750_MIMAT0004692 </td> <td> UP </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-452 </td> <td> MIMAT0001635 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
  <tr> <td> miR-17* </td> <td> MIMAT0000071 </td> <td> DOWN </td> <td align="right"> 0.015 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> </tr>
   </table>


## Checking results towards a curated list of microRNA´s

```r
	curatedfn = "curated_microRNA.csv"
  curatedtab  = read.table(paste("not_in_github", "/", curatedfn, sep=""), sep="\t", 
  												 header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                                     strip.white=TRUE, comment.char ="")
	curatedtab$miRNA = paste("hsa-", curatedtab$miRNA, sep="")
  # ad hoc fix
	curatedtab$miRNA = gsub("-mir-", "-miR-", curatedtab$miRNA)

	curatedtab$MIMAT = mimatmapping[curatedtab$miRNA]
	table(curatedtab$status)
```

```
## 
##   purple rejected    white   yellow 
##      302      969      484        4
```

```r
	curatedtab = curatedtab[!is.na(curatedtab$MIMAT),]
	table(curatedtab$status)
```

```
## 
##   purple rejected    white   yellow 
##      244      515      347        4
```


Adding the curated status to the results


```r
for(n in names(volinia_overlap))
{
	volinia_overlap[[n]]$curated = curatedtab[ match( rownames(volinia_overlap[[n]]), curatedtab$MIMAT ), "status"]
	voliniatab[[n]]$curated = curatedtab[ match( as.character(voliniatab[[n]]$MIMATID), curatedtab$MIMAT ), "status"]
}
x = order(volinia_overlap[["DCIS-normal"]]$FDR_merged)
print( xtable( volinia_overlap[["DCIS-normal"]][x,] , caption="With curated status", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = FALSE)
```

<table CELLPADDING=5>
<caption align="bottom"> With curated status </caption>
<tr> <th> miRNA_vol </th> <th> MIMATID_vol </th> <th> Direction_vol </th> <th> FDR_vol </th> <th> Direction_merged </th> <th> FDR_merged </th> <th> Direction_meta </th> <th> FDR_meta </th> <th> validated </th> <th> curated </th>  </tr>
  <tr> <td> miR-21 </td> <td> MIMAT0000076 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-497 </td> <td> MIMAT0002820 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-145 </td> <td> MIMAT0000437 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-193a-5p </td> <td> MIMAT0004614 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-378 </td> <td> MIMAT0000731_MIMAT0000732 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-155 </td> <td> MIMAT0000646 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-200c </td> <td> MIMAT0000617 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-425 </td> <td> MIMAT0001343_MIMAT0003393 </td> <td> UP </td> <td align="right"> 0.012 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-183 </td> <td> MIMAT0000261 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-200b </td> <td> MIMAT0000318 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-429 </td> <td> MIMAT0001536 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-96 </td> <td> MIMAT0000095 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-142-3p </td> <td> MIMAT0000434 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-143* </td> <td> MIMAT0004599 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-140-3p </td> <td> MIMAT0004597 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-142-5p </td> <td> MIMAT0000433 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-99a </td> <td> MIMAT0000097 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> let-7b </td> <td> MI0000063_MIMAT0000063 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> let-7c </td> <td> MI0000064_MIMAT0000064 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-193b </td> <td> MIMAT0002819 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-106b </td> <td> MIMAT0000680 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.004 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-145* </td> <td> MIMAT0004601 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td>  </td> </tr>
  <tr> <td> miR-342 </td> <td> MIMAT0000753 </td> <td> UP </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-100 </td> <td> MIMAT0000098 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.008 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-210 </td> <td> MIMAT0000267 </td> <td> DOWN </td> <td align="right"> 0.009 </td> <td> UP </td> <td align="right"> 0.022 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-22 </td> <td> MIMAT0000077 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.024 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-652 </td> <td> MIMAT0003322 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.041 </td> <td> DOWN </td> <td align="right"> 0.050 </td> <td> TRUE </td> <td> white </td> </tr>
  <tr> <td> miR-374b </td> <td> MIMAT0004955 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.051 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-19a </td> <td> MIMAT0000073 </td> <td> UP </td> <td align="right"> 0.008 </td> <td> UP </td> <td align="right"> 0.073 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-376c </td> <td> MIMAT0000720 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> DOWN </td> <td align="right"> 0.088 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-221 </td> <td> MIMAT0000278 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.112 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-361-5p </td> <td> MIMAT0000703 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.114 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-30d </td> <td> MIMAT0000245 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> DOWN </td> <td align="right"> 0.116 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-223 </td> <td> MIMAT0000280 </td> <td> UP </td> <td align="right"> 0.003 </td> <td> UP </td> <td align="right"> 0.124 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-185 </td> <td> MIMAT0000455 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td> UP </td> <td align="right"> 0.130 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-374a </td> <td> MIMAT0000727 </td> <td> UP </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.195 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-127-3p </td> <td> MIMAT0000446 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> DOWN </td> <td align="right"> 0.200 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-26b </td> <td> MIMAT0000083 </td> <td> UP </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.377 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-574-3p </td> <td> MIMAT0003239 </td> <td> DOWN </td> <td align="right"> 0.000 </td> <td> UP </td> <td align="right"> 0.378 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-20a </td> <td> MIMAT0000075 </td> <td> UP </td> <td align="right"> 0.007 </td> <td> UP </td> <td align="right"> 0.446 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-107 </td> <td> MIMAT0000104 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.468 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-125a </td> <td> MIMAT0000443 </td> <td> DOWN </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.782 </td> <td> UP </td> <td align="right"> 0.050 </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-423-5p </td> <td> MIMAT0004748 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td> UP </td> <td align="right"> 0.852 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td>  </td> </tr>
  <tr> <td> let-7d </td> <td> MI0000065_MIMAT0000065 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td> UP </td> <td align="right"> 0.887 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-15b </td> <td> MIMAT0000417 </td> <td> UP </td> <td align="right"> 0.001 </td> <td> DOWN </td> <td align="right"> 0.945 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-29c </td> <td> MIMAT0000681 </td> <td> UP </td> <td align="right"> 0.006 </td> <td> UP </td> <td align="right"> 0.951 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-15a </td> <td> MIMAT0000068 </td> <td> UP </td> <td align="right"> 0.013 </td> <td> DOWN </td> <td align="right"> 0.972 </td> <td>  </td> <td align="right">  </td> <td> FALSE </td> <td> white </td> </tr>
  <tr> <td> miR-182 </td> <td> MIMAT0000259 </td> <td> UP </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td> white </td> </tr>
  <tr> <td> miR-423-3p </td> <td> MIMAT0001340 </td> <td> DOWN </td> <td align="right"> 0.001 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td> white </td> </tr>
  <tr> <td> miR-106b* </td> <td> MIMAT0004672 </td> <td> DOWN </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td>  </td> </tr>
  <tr> <td> miR-32 </td> <td> MIMAT0000090 </td> <td> UP </td> <td align="right"> 0.002 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td> white </td> </tr>
  <tr> <td> miR-28-3p </td> <td> MIMAT0004502 </td> <td> DOWN </td> <td align="right"> 0.004 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td>  </td> </tr>
  <tr> <td> miR-340 </td> <td> MIMAT0000750_MIMAT0004692 </td> <td> UP </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td> white </td> </tr>
  <tr> <td> miR-452 </td> <td> MIMAT0001635 </td> <td> DOWN </td> <td align="right"> 0.005 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td> white </td> </tr>
  <tr> <td> miR-17* </td> <td> MIMAT0000071 </td> <td> DOWN </td> <td align="right"> 0.015 </td> <td>  </td> <td align="right">  </td> <td>  </td> <td align="right">  </td> <td>  </td> <td>  </td> </tr>
   </table>


The difference in curated status for the filtered out microRNAs and the ones kept as expressed, and the ones found as significantly differentially expressed between the "DCIS" and "normal" stat. AHUS only for the "expressed" and "nonexpressed". Both UCAM and AHUS in the significant. 



```r
#	table(curatedtab$status[  curatedtab$MIMAT %in% rownames(mergedtables[[1]]) ])

expressed =	table(curatedtab$status[curatedtab$MIMAT %in% filteredAHUSMIMAT])
nonexpressed=	table(curatedtab$status[curatedtab$MIMAT %in% unfilteredAHUSMIMAT[!unfilteredAHUSMIMAT %in% filteredAHUSMIMAT]])

genelist = intersect( rownames(mergedtables[["DCIS-normal"]])[mergedtables[["DCIS-normal"]]$adj.P.Val<0.05], 
					 rownames(metatables[["DCIS-normal"]]) )
significant = table(curatedtab$status[curatedtab$MIMAT %in% genelist])


print( xtable( rbind(nonexpressed,expressed,significant), caption="", digits=3), 
      comment = FALSE,
      type = "html",
      html.table.attributes="CELLPADDING=5",
      include.rownames = TRUE)
```

<table CELLPADDING=5>
<caption align="bottom">  </caption>
<tr> <th>  </th> <th> purple </th> <th> rejected </th> <th> white </th> <th> yellow </th>  </tr>
  <tr> <td align="right"> nonexpressed </td> <td align="right">   64 </td> <td align="right">   94 </td> <td align="right">  159 </td> <td align="right">    2 </td> </tr>
  <tr> <td align="right"> expressed </td> <td align="right">    3 </td> <td align="right">   24 </td> <td align="right">  145 </td> <td align="right">    1 </td> </tr>
  <tr> <td align="right"> significant </td> <td align="right">    1 </td> <td align="right">    5 </td> <td align="right">   46 </td> <td align="right">    1 </td> </tr>
   </table>


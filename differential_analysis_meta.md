
# Meta analysis on AHUS and UCAM miRNA microarrays datasets 
## Liliana Greger

## Data description
The  dataset  consists of :
* The AHUS datasets are Agilent_miRNA_v2_4470B consisting of  4 benign, 15 invasive, 18 normal, no DICS samples
* Agilent_miRNA_v3_4470C consists of  23 benign, 55 invasive, 70 normal and  8 DICS
* UCAM dataset consists of  8 benign, 1283 invasive, 116 normal and 10 DICS

AHUS outputs comes from Agilents Feature Extraction Software (AFE). UCAM data  has been received  as a preprocessed matrix as described in the Nature paper (Nature. 2013 May 16;497(7449):378-82). 

```r
library(AgiMicroRna)
library(Biobase)
library(MAMA)
```

```
## Loading required package: genefilter
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
## 
## Loading required package: metaMA
## Loading required package: SMVar
## 
## Attaching package: 'metaMA'
## 
## The following object is masked from 'package:genefilter':
## 
##     rowVars
## 
## Loading required package: multtest
## Loading required package: gtools
## 
## Attaching package: 'gtools'
## 
## The following object is masked from 'package:mgcv':
## 
##     scat
## 
## Loading required package: grid
## Loading required package: GeneMeta
## 
## Attaching package: 'MAMA'
## 
## The following objects are masked from 'package:GeneMeta':
## 
##     multExpFDR, zScoreFDR, zScorePermuted, zScores
```

```r
library(plyr)
#load("AHUS_agiMicroRNAdata.rdata")
if(!exists("inputisread"))
	source("read_input.r")
```

## Preprocessing -  summarization
Was done prior, see read_input.r.

```r
annot <- sampleannotation 
## use the UCAM data loaded earlier.
ucam = ucam_matrix
```
## Create the target files
We further annotate the samples using the annotation file.  No information on matched samples from UCAM.  

```r
## create ahus 3 target file 
ahus_v3_targets <- ahus_uRNAList$targets
ahus_3 <- as.matrix(ahus_v3_targets$FileName)
status <- c()
 for ( i in 1: length(ahus_v3_targets$FileName)) {
status = c(status, as.vector(annot[which(annot$datafile %in% ahus_3[i]), "tissue_type"]))}
var2 <- c()
for ( i in 1: length(ahus_v3_targets$FileName)) {
var2 = c(var2, as.vector(annot[which(annot$datafile %in% ahus_3[i]), "pam50_est"]))}
ihc <- c()
for ( i in 1: length(ahus_v3_targets$FileName)) {
ihc = c(ihc, as.vector(annot[which(annot$datafile %in% ahus_3[i]), "IHC"]))}
GErep <- gsub("benign", "1", status)
GErep <- gsub("DCIS", "2", GErep)
GErep <- gsub("invasive", "3", GErep)
GErep <- gsub("normal", "4", GErep)
ahus_v3_targets <- lapply(ahus_v3_targets, function(x) gsub("microRNA/AHUS/data/Agilent_miRNA_v3_4470C/", "", x))
ahus_v3_targets  <- sapply(ahus_v3_targets , function(x) sub("_.*", "", x))
ahus_v3_targets <- as.data.frame(ahus_v3_targets)
targets_ahus_v3 <- cbind(ahus_v3_targets ,status,  var2,ihc,  GErep) 
targets_ahus_v3$var2 = as.character(targets_ahus_v3$var2)
targets_ahus_v3$var2[which(targets_ahus_v3$status == "DCIS") ] = "DCIS" #set all DCIS to "DCIS" for the pam50 score since there are so few DCIS samples with pam50 score.
targets_ahus_v3$var2 = factor(targets_ahus_v3$var2)
targets_ahus_v3$ihc <- as.character.factor(targets_ahus_v3$ihc)
targets_ahus_v3$ihc[which(targets_ahus_v3$var2 == "DCIS") ] = "DCIS"
targets_ahus_v3$ihc <- as.factor(targets_ahus_v3$ihc)
##  create  ucam target file
status <- c()
for ( i in 1 : length(colnames(ucam))) {
status <- c(status, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$sample_id, value = "TRUE"), "tissue_type"]))}
GErep <- gsub("benign", "1", status)
GErep <- gsub("DCIS", "2", GErep)
GErep <- gsub("invasive", "3", GErep)
GErep <- gsub("normal", "4", GErep)
var2<- c()
for ( i in 1 : length(colnames(ucam))) {
var2<- c(var2, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$sample_id, value = "TRUE"), "pam50_est"]))}
ihc <- c()
for ( i in 1 : length(colnames(ucam))) {
ihc <- c(ihc, as.vector(annot[ annot$sample_id %in% grep( colnames(ucam)[i],annot$sample_id, value = "TRUE"), "IHC"]))}
targets_ucam = as.data.frame(cbind(FileName = colnames(ucam), status, var2,ihc,  GErep))
targets_ucam$var2 = as.character(targets_ucam$var2)
targets_ucam$var2[which(targets_ucam$status == "DCIS") ] = "DCIS" #set all DCIS to "DCIS" for the pam50 score since there are so few DCIS samples with pam50 score.
targets_ucam$var2 = factor(targets_ucam$var2)
targets_ucam$ihc <- as.character.factor(targets_ucam$ihc)
targets_ucam$ihc[which(targets_ucam$var2 == "DCIS") ] = "DCIS"
targets_ucam$ihc <- as.factor(targets_ucam$ihc)
```

## Creating the expression sets 
We then create the expression sets and later make sure the all datasets have common feature names in the exact order

```r
ahus_v3_eset =esetMicroRna(ahus_v3_filtered ,targets_ahus_v3 ,makePLOT=FALSE,verbose=TRUE)
```

```
## outPUT DATA: esetPROC 
## Features  Samples 
##      288      156 
## ------------------------------------------------------
```

```r
ucam_data <- targets_ucam
rownames(ucam_data) <- targets_ucam$FileName
pd <- new("AnnotatedDataFrame", data = ucam_data)
ucam_eset <- new("ExpressionSet", exprs = as.matrix(ucam) , phenoData = pd)
```
##  Conversion of   miRNA names to miRbase accession numbers
A microRNA to MIMATID mapping was made earlier (in read_input.r) . The UCAM data is already using MIMAT as ID.

```r
featureNames(ahus_v3_eset) = miRNA2MIMAT[match(featureNames(ahus_v3_eset), miRNA2MIMAT$AHUS), "MIMAT"]
```
## Find common features for both datasets

```r
ahus_v3_features <- featureNames(ahus_v3_eset) ## 266
ucam_features <- featureNames(ucam_eset)
all_features <- intersect(ahus_v3_features, ucam_features) ## 
ahus_v3_common_features <- ahus_v3_eset[featureNames(ahus_v3_eset)  %in% all_features , ]
ahus_v3_common_features <- ahus_v3_common_features[order( featureNames(ahus_v3_common_features)),]
ucam_common_features <- ucam_eset[featureNames(ucam_eset)  %in% all_features , ]
ucam_common_features <- ucam_common_features[order( featureNames(ucam_common_features)),]
```
## Perform meta analysis
We use the combined p value meta analysis using limma implemented in the MAMA package 
## First, prepare the ahus sub-dataseta for meta analysis

```r
ahus_v3_dcis_lumA <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "LumA")]
ahus_v3_dcis_lumA <- ahus_v3_dcis_lumA[ , ahus_v3_dcis_lumA$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_lumA))$var2 <- factor( pData(phenoData(ahus_v3_dcis_lumA))$var2) 
pData(phenoData(ahus_v3_dcis_lumA))$FileName <- factor( pData(phenoData(ahus_v3_dcis_lumA))$FileName)
pData(phenoData(ahus_v3_dcis_lumA))$status <- factor( pData(phenoData(ahus_v3_dcis_lumA))$status) 
ahus_v3_dcis_lumB <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "LumB")]
ahus_v3_dcis_lumB <- ahus_v3_dcis_lumB[ , ahus_v3_dcis_lumB$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_lumB))$var2 <- factor( pData(phenoData(ahus_v3_dcis_lumB))$var2) 
pData(phenoData(ahus_v3_dcis_lumB))$FileName <- factor( pData(phenoData(ahus_v3_dcis_lumB))$FileName) 
pData(phenoData(ahus_v3_dcis_lumB))$status <- factor( pData(phenoData(ahus_v3_dcis_lumB))$status) 
ahus_v3_dcis_Basal <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in% c("DCIS", "Basallike")]
ahus_v3_dcis_Basal<- ahus_v3_dcis_Basal[ , ahus_v3_dcis_Basal$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Basal))$var2 <- factor(pData(phenoData(ahus_v3_dcis_Basal))$var2) 
pData(phenoData(ahus_v3_dcis_Basal))$FileName <- factor(pData(phenoData(ahus_v3_dcis_Basal))$FileName)
pData(phenoData(ahus_v3_dcis_Basal))$status <- factor(pData(phenoData(ahus_v3_dcis_Basal))$status)
ahus_v3_dcis_Normal <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "Normallike")]
ahus_v3_dcis_Normal <- ahus_v3_dcis_Normal[ , ahus_v3_dcis_Normal$status %in% c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Normal))$var2 <- factor(pData(phenoData(ahus_v3_dcis_Normal))$var2) 
pData(phenoData(ahus_v3_dcis_Normal))$FileName <- factor(pData(phenoData(ahus_v3_dcis_Normal))$FileName) 
pData(phenoData(ahus_v3_dcis_Normal))$status <- factor(pData(phenoData(ahus_v3_dcis_Normal))$status) 
ahus_v3_dcis_Her2 <- ahus_v3_common_features[ , ahus_v3_common_features$var2 %in%  c("DCIS", "Her2")]
ahus_v3_normal_Her2<- ahus_v3_dcis_Her2[ , ahus_v3_dcis_Her2$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_dcis_Her2))$var2 <-factor( pData(phenoData(ahus_v3_dcis_Her2))$var2)  
pData(phenoData(ahus_v3_dcis_Her2))$FileName <- factor( pData(phenoData(ahus_v3_dcis_Her2))$FileName) 
pData(phenoData(ahus_v3_dcis_Her2))$status <- factor( pData(phenoData(ahus_v3_dcis_Her2))$status) 
ahus_v3_normal_DCIS <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("normal", "DCIS")]
pData(phenoData(ahus_v3_normal_DCIS))$status <- factor( pData(phenoData(ahus_v3_normal_DCIS))$status)# drop previous factor levels
pData(phenoData(ahus_v3_normal_DCIS))$FileName <- factor( pData(phenoData(ahus_v3_normal_DCIS))$FileName)# 
ahus_v3_normaltobenign <-  ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("normal", "benign")]
pData(phenoData(ahus_v3_normaltobenign ))$status <- factor( pData(phenoData(ahus_v3_normaltobenign ))$status)# drop previous factor levels
pData(phenoData(ahus_v3_normaltobenign ))$FileName <- factor( pData(phenoData(ahus_v3_normaltobenign ))$FileName)# 
ahus_v3_invasive_DCIS <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
pData(phenoData(ahus_v3_invasive_DCIS))$status <- factor( pData(phenoData(ahus_v3_invasive_DCIS))$status)# drop previous factor levels
pData(phenoData(ahus_v3_invasive_DCIS))$FileName <- factor( pData(phenoData(ahus_v3_invasive_DCIS))$FileName)# 
ahus3_dcis_HER2neg_ERneg_PGRneg  <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2neg_ERneg_PGRneg <- ahus3_dcis_HER2neg_ERneg_PGRneg[, ahus3_dcis_HER2neg_ERneg_PGRneg$ihc %in% c("DCIS", "HER2neg_ERneg_PGRneg")]
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$ihc) 
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$FileName) 
pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$status <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg))$status) 
ahus3_dcis_HER2neg_ERpos  <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2neg_ERpos <- ahus3_dcis_HER2neg_ERpos[, ahus3_dcis_HER2neg_ERpos$ihc %in% c("DCIS", "HER2neg_ERpos")]
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$ihc) 
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$FileName) 
pData(phenoData(ahus3_dcis_HER2neg_ERpos))$status <- factor( pData(phenoData(ahus3_dcis_HER2neg_ERpos))$status) 
ahus3_dcis_HER2pos_ERneg   <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2pos_ERneg  <- ahus3_dcis_HER2pos_ERneg [, ahus3_dcis_HER2pos_ERneg $ihc %in% c("DCIS", "HER2pos_ERneg")]
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$ihc) 
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$FileName) 
pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$status <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERneg ))$status) 
ahus3_dcis_HER2pos_ERpos   <- ahus_v3_common_features[ , ahus_v3_common_features$status %in%  c("invasive", "DCIS")]
ahus3_dcis_HER2pos_ERpos  <- ahus3_dcis_HER2pos_ERpos [, ahus3_dcis_HER2pos_ERpos $ihc %in% c("DCIS", "HER2pos_ERpos")]
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$ihc <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$ihc) 
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$FileName <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$FileName) 
pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$status <- factor( pData(phenoData(ahus3_dcis_HER2pos_ERpos ))$status)
```
## Prepare ucam sub-datasets for meta analysis


```r
ucam_normalDCIS<- ucam_common_features[, ucam_common_features$status %in% c("normal", "DCIS") ]   
pData(phenoData(ucam_normalDCIS))$status  <- factor(pData(phenoData(ucam_normalDCIS))$status )
pData(phenoData(ucam_normalDCIS))$FileName  <- factor(pData(phenoData(ucam_normalDCIS))$FileName )
ucam_invasiveDCIS<- ucam_common_features[, ucam_common_features$status %in% c("invasive", "DCIS") ]
pData(phenoData(ucam_invasiveDCIS))$status  <- factor(pData(phenoData(ucam_invasiveDCIS))$status )               
pData(phenoData(ucam_invasiveDCIS))$FileName  <- factor(pData(phenoData(ucam_invasiveDCIS))$FileName )  
ucam_normalTobenign<- ucam_common_features[, ucam_common_features$status %in% c("benign", "normal") ]
pData(phenoData(ucam_normalTobenign))$status  <- factor(pData(phenoData(ucam_normalTobenign))$status )               
pData(phenoData(ucam_normalTobenign))$FileName  <- factor(pData(phenoData(ucam_normalTobenign))$FileName )
ucam_dcis_lumA <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "LumA") ]
ucam_dcis_lumA <- ucam_dcis_lumA[ , ucam_dcis_lumA$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_lumA))$var2  <- factor(pData(phenoData(ucam_dcis_lumA))$var2 )               
pData(phenoData(ucam_dcis_lumA))$FileName  <- factor(pData(phenoData(ucam_dcis_lumA))$FileName )
pData(phenoData(ucam_dcis_lumA))$status <- factor(pData(phenoData(ucam_dcis_lumA))$status)
ucam_dcis_lumB <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "LumB") ]
ucam_dcis_lumB <- ucam_dcis_lumB[ , ucam_dcis_lumB$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_lumB))$var2  <- factor(pData(phenoData(ucam_dcis_lumB))$var2 )               
pData(phenoData(ucam_dcis_lumB))$FileName  <- factor(pData(phenoData(ucam_dcis_lumB))$FileName )
pData(phenoData(ucam_dcis_lumB))$status <- factor(pData(phenoData(ucam_dcis_lumB))$status)
ucam_dcis_Her2 <- ucam_common_features[, ucam_common_features$var2 %in% c("DCIS", "Her2") ]
ucam_dcis_Her2 <- ucam_dcis_Her2[ , ucam_dcis_Her2$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Her2))$var2  <- factor(pData(phenoData(ucam_dcis_Her2))$var2 )               
pData(phenoData(ucam_dcis_Her2))$FileName  <- factor(pData(phenoData(ucam_dcis_Her2))$FileName )
pData(phenoData(ucam_dcis_Her2))$status  <- factor(pData(phenoData(ucam_dcis_Her2))$status )
ucam_dcis_Basal <- ucam_common_features[ , ucam_common_features$var2 %in% c("DCIS", "Basallike")]
ucam_dcis_Basal<- ucam_dcis_Basal[ , ucam_dcis_Basal$status %in%  c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Basal))$var2 <- factor(pData(phenoData(ucam_dcis_Basal))$var2) 
pData(phenoData(ucam_dcis_Basal))$FileName <- factor(pData(phenoData(ucam_dcis_Basal))$FileName)
pData(phenoData(ucam_dcis_Basal))$status <- factor(pData(phenoData(ucam_dcis_Basal))$status)
ucam_dcis_Normal <- ucam_common_features[ , ucam_common_features$var2 %in%  c("DCIS", "Normallike")]
ucam_dcis_Normal <- ucam_dcis_Normal[ , ucam_dcis_Normal$status %in% c("invasive", "DCIS")]
pData(phenoData(ucam_dcis_Normal))$var2 <- factor(pData(phenoData(ucam_dcis_Normal))$var2) 
pData(phenoData(ucam_dcis_Normal))$FileName <- factor(pData(phenoData(ucam_dcis_Normal))$FileName) 
pData(phenoData(ucam_dcis_Normal))$status <- factor(pData(phenoData(ucam_dcis_Normal))$status) 
ucam_dcis_HER2neg_ERneg_PGRneg  <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2neg_ERneg_PGRneg <- ucam_dcis_HER2neg_ERneg_PGRneg[, ucam_dcis_HER2neg_ERneg_PGRneg$ihc %in% c("DCIS", "HER2neg_ERneg_PGRneg")]
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$ihc <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$ihc) 
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$FileName <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$FileName) 
pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$status <- factor( pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg))$status) 
ucam_dcis_HER2neg_ERpos  <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2neg_ERpos <- ucam_dcis_HER2neg_ERpos[, ucam_dcis_HER2neg_ERpos$ihc %in% c("DCIS", "HER2neg_ERpos")]
pData(phenoData(ucam_dcis_HER2neg_ERpos))$ihc <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$ihc) 
pData(phenoData(ucam_dcis_HER2neg_ERpos))$FileName <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$FileName) 
pData(phenoData(ucam_dcis_HER2neg_ERpos))$status <- factor( pData(phenoData(ucam_dcis_HER2neg_ERpos))$status) 
ucam_dcis_HER2pos_ERneg   <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2pos_ERneg  <- ucam_dcis_HER2pos_ERneg [, ucam_dcis_HER2pos_ERneg $ihc %in% c("DCIS", "HER2pos_ERneg")]
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$ihc <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$ihc) 
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$FileName <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$FileName) 
pData(phenoData(ucam_dcis_HER2pos_ERneg ))$status <- factor( pData(phenoData(ucam_dcis_HER2pos_ERneg ))$status) 
ucam_dcis_HER2pos_ERpos   <- ucam_common_features[ , ucam_common_features$status %in%  c("invasive", "DCIS")]
ucam_dcis_HER2pos_ERpos  <- ucam_dcis_HER2pos_ERpos [, ucam_dcis_HER2pos_ERpos $ihc %in% c("DCIS", "HER2pos_ERpos")]
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$ihc <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$ihc) 
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$FileName <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$FileName) 
pData(phenoData(ucam_dcis_HER2pos_ERpos ))$status <- factor( pData(phenoData(ucam_dcis_HER2pos_ERpos ))$status)
```

## Prepare the combined sub-datatests


```r
all_normal_DCIS <- new("MetaArray", GEDM = list( exprs(ahus_v3_normal_DCIS), exprs(ucam_normalDCIS)), clinical = list(pData(phenoData(ahus_v3_normal_DCIS)), pData(phenoData(ucam_normalDCIS) )) , datanames = c( "ahus_v3_4470C_normal_DCIS", "ucam_normal_DCIS"))
all_invasive_DCIS <- new("MetaArray", GEDM = list( exprs(ahus_v3_invasive_DCIS), exprs(ucam_invasiveDCIS )), clinical = list(pData(phenoData(ahus_v3_invasive_DCIS)), pData(phenoData(ucam_invasiveDCIS ) )) , datanames = c( "ahus_v3_4470C_invasive_DCIS", "ucam_invasive_DCIS"))
all_normnal_to_benign <- new("MetaArray", GEDM = list( exprs(ahus_v3_normaltobenign), exprs(ucam_normalTobenign )), clinical = list(pData(phenoData(ahus_v3_normaltobenign)), pData(phenoData(ucam_normalTobenign ) )) , datanames = c( "ahus_v3_4470C_invasive_DCIS", "ucam_normalTobenign"))
all_dcis_LumA  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_lumA), exprs(ucam_dcis_lumA )), clinical = list(pData(phenoData(ahus_v3_dcis_lumA)), pData(phenoData(ucam_dcis_lumA ) )) , datanames = c( "ahus_v3_dcis_lumA", "ucam_dcis_lumA"))
all_dcis_LumB  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_lumB), exprs(ucam_dcis_lumB )), clinical = list(pData(phenoData(ahus_v3_dcis_lumB)), pData(phenoData(ucam_dcis_lumB ) )) , datanames = c( "ahus_v3_dcis_lumB", "ucam_dcis_lumB"))
all_dcis_Her2  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Her2), exprs(ucam_dcis_Her2 )), clinical = list(pData(phenoData(ahus_v3_dcis_Her2)), pData(phenoData(ucam_dcis_Her2 ) )) , datanames = c( "ahus_v3_dcis_lumB", "ucam_dcis_Her2"))
all_dcis_Normal  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Normal),exprs(ucam_dcis_Normal )), clinical =list(pData(phenoData(ahus_v3_dcis_Normal)), pData(phenoData(ucam_dcis_Normal ))) , datanames = c( "ahus_v3_dcis_Normal", "ucam_dcis_Normal"))
all_dcis_Basal  <- new("MetaArray", GEDM = list( exprs(ahus_v3_dcis_Basal),exprs(ucam_dcis_Basal )), clinical = list(pData(phenoData(ahus_v3_dcis_Basal)),pData(phenoData(ucam_dcis_Basal) )) , datanames = c( "ahus_v3_dcis_Basal", "ucam_dcis_Basal"))
all_dcis_HER2neg_ERneg_PGRneg <- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2neg_ERneg_PGRneg),exprs(ucam_dcis_HER2neg_ERneg_PGRneg )), clinical =list(pData(phenoData(ahus3_dcis_HER2neg_ERneg_PGRneg)), pData(phenoData(ucam_dcis_HER2neg_ERneg_PGRneg ))) , datanames = c( "ahus3_dcis_HER2neg_ERneg_PGRneg", "ucam_dcis_HER2neg_ERneg_PGRneg"))
all_dcis_HER2pos_ERpos<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2pos_ERpos),exprs(ucam_dcis_HER2pos_ERpos )), clinical = list(pData(phenoData(ahus3_dcis_HER2pos_ERpos)), pData(phenoData(ucam_dcis_HER2pos_ERpos))) , datanames = c( "ahus3_dcis_HER2pos_ERpos", "ucam_dcis_HER2pos_ERpos"))
all_dcis_HER2pos_ERneg<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2pos_ERneg),exprs(ucam_dcis_HER2pos_ERneg )), clinical = list(pData(phenoData(ahus3_dcis_HER2pos_ERneg)), pData(phenoData(ucam_dcis_HER2pos_ERneg))) , datanames = c( "ahus3_dcis_HER2pos_ERneg", "ucam_dcis_HER2pos_ERneg"))
all_dcis_HER2neg_ERpos<- new("MetaArray", GEDM = list(exprs(ahus3_dcis_HER2neg_ERpos),exprs(ucam_dcis_HER2neg_ERpos )), clinical = list(pData(phenoData(ahus3_dcis_HER2neg_ERpos)), pData(phenoData(ucam_dcis_HER2neg_ERpos))) , datanames = c( "ahus3_dcis_HER2neg_ERpos", "ucam_dcis_HER2neg_ERpos"))
```
## Meta analysis using limma


```r
outputdir = "output-analysis"
if(!file.exists(outputdir))
	dir.create(outputdir)
difflistdir = paste(outputdir, "/difflists_meta", sep="")
if(!file.exists(difflistdir))
		dir.create(difflistdir)

results_pval1 <- metaMA(all_normal_DCIS, varname = "status", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/DCIS_to_normal_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_invasive_DCIS, varname = "status", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/DCIS_to_invasive_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_normnal_to_benign, varname = "status", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/benign_to_normal_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_dcis_LumA, varname = "var2", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_lumA_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_dcis_LumB, varname = "var2", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_lumB_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_dcis_Her2, varname = "var2", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_her2_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_dcis_Normal, varname = "var2", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_normallike_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(all_dcis_Basal, varname = "var2", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/basallike_to_dcis_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(  all_dcis_HER2neg_ERneg_PGRneg, varname = "ihc", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_HER2neg_ERneg_PGRneg_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(  all_dcis_HER2neg_ERpos, varname = "ihc", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_HER2neg_ERpos_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(  all_dcis_HER2pos_ERneg, varname = "ihc", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_HER2pos_ERneg_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
results_pval1 <- metaMA(  all_dcis_HER2pos_ERpos, varname = "ihc", moderated = "limma",   which = "pval", BHth = 1) 
```

```
##   DE  IDD Loss  IDR  IRR 
##  265    0    0    0    0
```

```r
test <- data.frame(MIMAT = c(), miRNAname = c(),FC=c())
for (i in c(results_pval1$Meta) ){
DE <- results_pval1$gene.names[i]
feature <-  miRNA2MIMAT[DE, "preferredname"]
test_st <-  results_pval1$TestStatistic[i]
pval <- 2*(1-pnorm(abs(test_st)))
test2 <- data.frame(MIMAT = DE,  miRNAname = paste(feature, sep = ",", collapse = ","), TestStatistic = test_st, pval)
test <- rbind(test, test2)
}
test$adj.P.Val=p.adjust(test$pval, method="BH")
test = test[order(test$pval),]
write.table(test, "output-analysis/difflists_meta/dcis_to_HER2pos_ERpos_limma.txt", sep = "\t", quote = FALSE, row.names=FALSE)
```










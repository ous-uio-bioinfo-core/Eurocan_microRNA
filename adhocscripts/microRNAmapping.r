

# ad hoc script to make a master miRNA to MIMAT mapping for the combined data set, UCAM and AHUS.
# Annotations from the Agilent feature extraction file, mirBASE(alias.txt) and a few manual lookups are used.

# ahus_MIMAT_miRNA_key.txt
# UCAM_genes.txt
# alias.txt

# annotation from mirBASE.
tmpannotdir = "./"
annotdir = "../input"
aliases <- read.delim(paste(tmpannotdir, "/aliases.txt", sep=""), header = FALSE, sep = "") ## file from miRbase
aliases = aliases[grepl("MIMAT", aliases$V1), ]# only use MIMAT id and not MI
library(data.table)
geneannot <- data.table(aliases, key="V1")
geneannot <- geneannot[, list(V2 = unlist(strsplit(as.character(V2), ";"))), by=V1]
geneannot <- as.matrix(geneannot)
geneannot <- as.data.frame(geneannot)
names(geneannot) = c("MIMAT", "miRNA")
geneannot[,"miRNA"] = as.character(geneannot[,"miRNA"])
geneannot[,"MIMAT"] = as.character(geneannot[,"MIMAT"])
geneannot = geneannot[!duplicated(geneannot),]


# Make our master miRNA to MIMAT mapping.
miRNA2MIMAT = read.table(paste(tmpannotdir, "/", "ahus_MIMAT_miRNA_key.txt", sep=""), header = TRUE, sep="\t", 
												stringsAsFactors=FALSE)
names(miRNA2MIMAT)[1]="AHUS"
rownames(miRNA2MIMAT)=miRNA2MIMAT$MIMAT

### Map UCAM miRNA #######

ucam_genes = read.table(paste(tmpannotdir, "/", "UCAM_genes.txt", sep=""), header = FALSE, sep="\t", 
												stringsAsFactors=FALSE)[,1]
miRNA2MIMAT$UCAM=NA
a = miRNA2MIMAT$AHUS %in% ucam_genes
miRNA2MIMAT$UCAM[a] = miRNA2MIMAT$AHUS[a]
remainder = ucam_genes[!ucam_genes %in% miRNA2MIMAT$UCAM]

# try to map the reminder ucam using alias.txt, and add those that match to one MIMAT to master.
ucamtab = geneannot[geneannot$miRNA %in% remainder, ]
ucamtab = ucamtab[!duplicated(ucamtab), ]
tmp = table(ucamtab$miRNA)
ucamtab$miRNAmatch = tmp[match(ucamtab$miRNA, names(tmp))]
tmp = ucamtab[ ucamtab$miRNAmatch==1 & ucamtab$MIMAT %in% miRNA2MIMAT$MIMAT, c("miRNA", "MIMAT")]
miRNA2MIMAT[tmp$MIMAT, "UCAM"] = tmp$miRNA
tmp = ucamtab[ucamtab$miRNAmatch==1 & !(ucamtab$MIMAT %in% miRNA2MIMAT$MIMAT), c("miRNA", "MIMAT")]
miRNA2MIMAT = rbind(miRNA2MIMAT, data.frame(MIMAT=tmp$MIMAT, UCAM=tmp$miRNA, AHUS=NA, row.names=tmp$MIMAT))

#ucamtab[ucamtab$miRNAmatch!=1,]

# ad-hoc mmaping done on the above few
con <- textConnection(
"miRNA	MIMAT
hsa-miR-500	MIMAT0004773
hsa-miR-550	MIMAT0004800
kshv-miR-K12-12	MIMAT0003712
hsa-miR-3190-3p	MIMAT0022839")
#hsa-miR-550*	MIMAT0003257
#hsa-miR-500*	MIMAT0002871
adhocmap = read.table(con, header = TRUE, sep="\t", 
												 stringsAsFactors=FALSE)
tmp = adhocmap[ adhocmap$MIMAT %in% miRNA2MIMAT$MIMAT, c("miRNA", "MIMAT")]
miRNA2MIMAT[tmp$MIMAT, "UCAM"] = tmp$miRNA
tmp = adhocmap[!(adhocmap$MIMAT %in% miRNA2MIMAT$MIMAT), c("miRNA", "MIMAT")]
miRNA2MIMAT = rbind(miRNA2MIMAT, data.frame(MIMAT=tmp$MIMAT, UCAM=tmp$miRNA, AHUS=NA, row.names=tmp$MIMAT))


####### Map Volinia miRNA

genes1 = read.table(file=paste(annotdir, "/", "volinia-invasive-vs-DCIS.txt", sep=""),
																					 sep="\t", check.names=TRUE, 
																					 stringsAsFactors=FALSE, header=TRUE)[,1]
genes2 = read.table(file=paste(annotdir, "/", "volinia-DCIS-vs-normal.txt", sep=""),
																					sep="\t", check.names=TRUE, 
																					stringsAsFactors=FALSE, header=TRUE)[,1]
volinia_genes = c(genes1,genes2)
names(volinia_genes) = paste("hsa-", volinia_genes, sep="")
#volinia_genes_hsa = paste("hsa-", volinia_genes, sep="")

miRNA2MIMAT$volinia=NA
a = miRNA2MIMAT$AHUS %in% names(volinia_genes)
miRNA2MIMAT$volinia[a] = volinia_genes[miRNA2MIMAT$AHUS[a]]
remainder = volinia_genes[!volinia_genes %in% miRNA2MIMAT$volinia]

# try to map the reminder ucam using alias.txt, and add those that match to one MIMAT to master.
voliniatab = geneannot[geneannot$miRNA %in% names(remainder), ]
voliniatab = voliniatab[!duplicated(voliniatab), ]
tmp = voliniatab[ voliniatab$MIMAT %in% miRNA2MIMAT$MIMAT, c("miRNA", "MIMAT")]
miRNA2MIMAT[tmp$MIMAT, "volinia"] = remainder[tmp$miRNA]
tmp = voliniatab[!(voliniatab$MIMAT %in% miRNA2MIMAT$MIMAT), c("miRNA", "MIMAT")] # no matches
remainder = volinia_genes[!volinia_genes %in% miRNA2MIMAT$volinia]
#write.table(remainder, file="unmapped_volinia.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# ad-hoc mmaping done on the some of the above few
con <- textConnection(
	"miRNA	MIMAT
miR-16(2)	MIMAT0000069
miR-19b(2)	MIMAT0000074
miR-26a(2)	MIMAT0000082
miR-92a(2)	MIMAT0000092
miR-29b(2)	MIMAT0000100
miR-181a(2)	MIMAT0000256
miR-218(2)	MIMAT0000275
miR-125b(2)	MIMAT0000423
miR-376a-3p(2)	MIMAT0000729
miR-324	MIMAT0000761")
adhocmap = read.table(con, header = TRUE, sep="\t", 
											stringsAsFactors=FALSE)
tmp = adhocmap[ adhocmap$MIMAT %in% miRNA2MIMAT$MIMAT, c("miRNA", "MIMAT")]
miRNA2MIMAT[tmp$MIMAT, "volinia"] = tmp$miRNA
tmp = adhocmap[!(adhocmap$MIMAT %in% miRNA2MIMAT$MIMAT), c("miRNA", "MIMAT")] # none


# based on the mirBASE, attach the miRNA name with arm info, -3p or -5p
miRNA2MIMAT$arm=NA
miRNA2MIMAT$armcount=NA
miRNA2MIMAT$allarms = NA
for(i in 1:nrow(miRNA2MIMAT))
{
	id = miRNA2MIMAT$MIMAT[i]
	a = geneannot$MIMAT %in% id
	tmpnames = geneannot[a, "miRNA"]
	b = grepl("-5p|-3p" , tmpnames)
	#miRNA2MIMAT$armcount[i] = sum(b)
	#miRNA2MIMAT$allarms[i] = paste(tmpnames[b], collapse=", ")
	if(sum(b)==1)
	{
		miRNA2MIMAT[i, "arm"] = tmpnames[b]
	}else if(sum(b)>1){ # if MIMAT  matches to many arms, check if one of them are given and is the same in both data set.
		if(miRNA2MIMAT$UCAM[i]==miRNA2MIMAT$AHUS[i] & miRNA2MIMAT$AHUS[i] %in%  tmpnames[b])
			miRNA2MIMAT[i, "arm"] =  miRNA2MIMAT$AHUS[i]
	}
}

#miRNA2MIMAT[miRNA2MIMAT$armcount>1,]


# nice order
miRNA2MIMAT = miRNA2MIMAT[, c("MIMAT", "arm", "AHUS", "UCAM", "volinia")]
# This should map miRNA to exactly one MIMAT.
write.table(miRNA2MIMAT, file=paste(annotdir, "/miRNA2MIMAT.txt", sep=""), 
						row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

a=rowSums(is.na(miRNA2MIMAT[,c("AHUS", "UCAM")]))==0
table(miRNA2MIMAT$AHUS[a]==miRNA2MIMAT$UCAM[a]) # miRNA that were rescude based on this mapping compared to a straight miRNA name mapping.
miRNA2MIMAT[a & is.na(miRNA2MIMAT$arm),]













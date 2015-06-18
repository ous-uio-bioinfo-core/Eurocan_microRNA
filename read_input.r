


# Basic script to read in the input needed for the analyses performed in the "Volinia check article"
#
# Two different data sets of microRNA Agilent microarray data is used.
# UCAM data 
# AHUS data
#
# Sampleannotation from both sources was combined in one file.
#
# The microRNA names are mapped to MIMATIDs in order to combine results.

library(plyr)
annotdir = "input"
datadir = "not_in_github"
#library(AgiMicroRna)

AHUS_binary_data_file="AHUS_agiMicroRNAdata.rdata"
UCAM_data_file = "Agilent_ncRNA_60k_normalised_miRNA_expression_ORIGINAL.txt"
sampleannotation_file="sampleannotation_validatingvolinia.txt"
mimatmapping_file = "MIMATID_conversion.txt"

# the miRNA to MIMAT mapping with preferred arm name if availible, else use AHUS version
miRNA2MIMAT = read.table(file=paste(annotdir, "/miRNA2MIMAT.txt", sep=""), 	header=TRUE ,sep="\t", stringsAsFactors=FALSE)
rownames(miRNA2MIMAT)=miRNA2MIMAT$MIMAT
miRNA2MIMAT$preferredname=miRNA2MIMAT$arm
miRNA2MIMAT$preferredname[is.na(miRNA2MIMAT$preferredname)]=miRNA2MIMAT$AHUS[is.na(miRNA2MIMAT$preferredname)]

# read sampleannotation table
sampleannotation=read.table(paste(annotdir,  "/", sampleannotation_file, sep=""),
                            sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                            strip.white=TRUE, comment.char ="")
row.names(sampleannotation) = sampleannotation$sample_id
sampleannotation$IHC[sampleannotation$IHC==""] = "unknown"

## Preprocessing -  summarization
#First we load  and preprocess the AHUS data. The use only one of the AHUS dataset which has DCIS samples. 
#We choose to use  RMA summarisation of  the signal without background correction as described in BMC Research Notes 2010, 3:18.
#We also read the processed and already normalized.
#W# load the AHUS data as a uRNAlist object. Prepared earlier from Agilent Feature Extraction Result Files
# We then filter out the ahus data from control spots,  miRNA whose gMeanSignal is close to the 
# expression of the negative controls(25\% limit) and  miRNA which are not expressed. 
# We set up a limit of 75\% of the miRNA which must remain in at least one experimental condition. 
# 288 features for ahus_v3 dataset are left  out of 961.
load(paste(datadir, "/AHUS_rma_uRNAList.rdata", sep=""))
a =match(ahus_uRNAList$targets$FileName, sampleannotation$datafile)
filtertargets = sampleannotation[a , c("sample_id", "tissue_type")]
names(filtertargets) = c("FileName", "Treatment")
filtertargets$GErep = as.numeric(factor(filtertargets$Treatment))
ahus_v3_filtered =filterMicroRna(ahus_uRNAList , AHUS_agiMicroRNAdata$Agilent_miRNA_v3_4470C, control = TRUE, 
                                 IsGeneDetected=TRUE,limIsGeneDetected=75,limNEG=25, 
                                 targets =filtertargets , 
                                 verbose=TRUE,writeout=FALSE)
ahus_matrix = exprs(esetMicroRna(ahus_v3_filtered, filtertargets))
colnames(ahus_matrix) = filtertargets$FileName
ahus_matrix = ahus_matrix[rownames(ahus_matrix) %in% miRNA2MIMAT$AHUS, ]#mimatmapping
rownames(ahus_matrix) = miRNA2MIMAT[match(rownames(ahus_matrix),miRNA2MIMAT$AHUS), "MIMAT"]


# read the UCAM data, using the providers version
ucam_matrix = read.table(paste(datadir, "/", UCAM_data_file, sep=""), header = TRUE, sep="\t", stringsAsFactors=FALSE)
ucam_matrix = ucam_matrix[ucam_matrix[,2] %in% miRNA2MIMAT$UCAM, ]

# UCAM data has been already filtered - 823 features.
# But some microRNA are measured with several different probes, use only one probe per microRNA,
# chose the one with highest sd across samples.
rowsd = apply(ucam_matrix[,-c(1,2)], MARGIN=1, FUN=sd)
ucam_matrix = ucam_matrix[order(rowsd, decreasing=TRUE), ]
ucam_matrix = ucam_matrix[!duplicated(ucam_matrix[,2]), ]


ucam_microrna_names_filtered=ucam_matrix[,2]
rownames(ucam_matrix) =  miRNA2MIMAT[match(ucam_matrix[,2],miRNA2MIMAT$UCAM), "MIMAT"]
ucam_matrix = as.matrix(ucam_matrix[, -c(1,2)])
colnames(ucam_matrix) = paste("UCAM_", colnames(ucam_matrix), sep="") # another naming scheme is used in the annotation used here.

# A few of the samples are without sample annotation, probably some kind of control samples. We filter them out.
ucam_matrix = ucam_matrix[, colnames(ucam_matrix) %in% sampleannotation$sample_id]



# ##########klipp
# # Make a vecor of all common MIMAT before filtering.
# ucam_genes = read.table(paste(annotdir, "/", "UCAM_genes.txt", sep=""), header = FALSE, sep="\t", 
# 												stringsAsFactors=FALSE)[,1] # The IDs of all genes from the UCAM platform, also the ones filtered out in the data matrix.
# unfilteredUCAMMIMAT =  mimatmapping[ucam_genes]
# #rm(ucam_genes)
# unfilteredUCAMMIMAT = unfilteredUCAMMIMAT[!is.na(unfilteredUCAMMIMAT)]
# unfilteredUCAMMIMAT = unfilteredUCAMMIMAT[!duplicated(unfilteredUCAMMIMAT)]
# filteredUCAMMIMAT =   mimatmapping[ucam_microrna_names_filtered]
# filteredUCAMMIMAT = filteredUCAMMIMAT[!is.na(filteredUCAMMIMAT)]
# 
# filteredAHUSMIMAT =  rownames(ahus_matrix)
# unfilteredAHUSMIMAT =  mimatmapping[ahus_uRNAList$genes$GeneName]
# unfilteredAHUSMIMAT = unfilteredAHUSMIMAT[!is.na(unfilteredAHUSMIMAT)]
# unfilteredcommonMIMAT = intersect(unfilteredAHUSMIMAT, unfilteredUCAMMIMAT)
# 
# # # a back to name mapping for the two data sets
# 
# 
# # make a mapping back to the microRNA name. But use the arm name when availible.
# unfilteredcommonname = names(unfilteredAHUSMIMAT[unfilteredAHUSMIMAT %in% unfilteredcommonMIMAT])
# names(unfilteredcommonname) = mimatmapping[unfilteredcommonname]
# 
# for(i in 1:length(unfilteredcommonname))
# {
# 	id = names(unfilteredcommonname)[i]
# 	a = geneannot$ID %in% id
# 	tmpnames = geneannot[a, 1]
# 	b = grepl("-5p|-3p" , tmpnames)
# 	if(sum(b)==1)
# 		unfilteredcommonname[i] = tmpnames[b]	
# }
# ########################klipp

# filter data based on common MIMATS,
filteredcommonMIMAT = intersect(rownames(ucam_matrix), rownames(ahus_matrix))
ucam_matrix = ucam_matrix[filteredcommonMIMAT, ]
ahus_matrix = ahus_matrix[filteredcommonMIMAT, ]
common_matrix= cbind(ahus_matrix, ucam_matrix)

inputisread = TRUE

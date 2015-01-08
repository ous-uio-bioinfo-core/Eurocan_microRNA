


# Basic script to read in the input needed for the analyses performed in the "Volinia check article"
#
# Two different data sets of microRNA Agilent microarray data is used.
# UCAM data 
# AHUS data
#
# Sampleannotation from both sources was combined in one file.
#
# The microRNA names are mapped to MIMATIDs in order to combine results.


annotdir = "input"
datadir = "not_in_github"
#library(AgiMicroRna)

AHUS_binary_data_file="AHUS_agiMicroRNAdata.rdata"
UCAM_data_file = "Agilent_ncRNA_60k_normalised_miRNA_expression_ORIGINAL.txt"
sampleannotation_file="sampleannotation_validatingvolinia.txt"
mimatmapping_file = "MIMATID_conversion.txt"


# make MIMAT mapping
mimattab= read.table(file=paste(annotdir, "/", mimatmapping_file, sep=""),
                     sep="\t", check.names=TRUE, comment.char="",
                     stringsAsFactors=FALSE, header=TRUE)
mimatmapping = mimattab[,2]
names(mimatmapping)= mimattab[,1]
rm(mimattab)

## make mapping from micorRNA name to MIMAT id, when several MIMAT matches, only the first alphabetically will be used.
# aliases <- read.delim(paste(annotdir ,"/aliases.txt", sep=""), header = FALSE, sep = "") ## file from miRbase
# colnames(aliases) <- c("MIMAT", "name")
# library(data.table)
# df <- data.table(aliases, key="MIMAT")
# df <- df[, list(name = unlist(strsplit(as.character(name), ";"))), by=MIMAT]
# df = as.data.frame(as.matrix(df), stringsAsFactors=FALSE)
# df = df[ grepl("MIMAT", df$MIMAT), ]
# df = df[order(df$name, df$MIMAT),]
# df = df[!duplicated(df$name),]
# mimatmapping = df$MIMAT
# names(mimatmapping)= df$name

# read sampleannotation table
sampleannotation=read.table(paste(annotdir,  "/", sampleannotation_file, sep=""),
                            sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
                            strip.white=TRUE, comment.char ="")
row.names(sampleannotation) = sampleannotation$sample_id
sampleannotation$IHC[sampleannotation$IHC==""] = "unknown"
#load(paste(datadir, "/", AHUS_binary_data_file, sep=""))
#ahus_uRNAList  = rmaMicroRna(AHUS_agiMicroRNAdata$Agilent_miRNA_v3_4470C, normalize = TRUE, background = FALSE)
#save(ahus_uRNAList, file="inputdata/AHUS_rma_uRNAList.rdata")

# load the AHUS data as a uRNAlist object. Prepared earlier from Agilent Feature Extraction Result Files
load(paste(datadir, "/AHUS_rma_uRNAList.rdata", sep=""))
# Format and filter to take out weak/unreliable/unmappable probes
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
ahus_matrix = ahus_matrix[rownames(ahus_matrix) %in% names(mimatmapping), ]#mimatmapping
rownames(ahus_matrix) = mimatmapping[rownames(ahus_matrix)]


# read the UCAM data, using the providers version
ucam_matrix = read.table(paste(datadir, "/", UCAM_data_file, sep=""), header = TRUE, sep="\t", stringsAsFactors=FALSE)
ucam_matrix = ucam_matrix[ucam_matrix[,2] %in%  names(mimatmapping),]

# some microRNA are measured with several different probes, use only one probe per microRNA, chose the one with highest sd across samples.
rowsd = apply(ucam_matrix[,-c(1,2)], MARGIN=1, FUN=sd)
ucam_matrix = ucam_matrix[order(rowsd, decreasing=TRUE), ]
ucam_matrix = ucam_matrix[!duplicated(ucam_matrix[,2]), ]

# debug
# x =  mimatmapping[ucam_matrix[,2]]
# y = duplicated(x) | duplicated(x, fromLast = TRUE)
# ucam_matrix[y,2]

rownames(ucam_matrix) = mimatmapping[ucam_matrix[,2]]
ucam_matrix = as.matrix(ucam_matrix[, -c(1,2)])
colnames(ucam_matrix) = paste("UCAM_", colnames(ucam_matrix), sep="") # another naming scheme is used in the annotation I made.

# A few of the samples are without sample annotation, probably some kind of control samples. I filter them out.
ucam_matrix = ucam_matrix[, colnames(ucam_matrix) %in% sampleannotation$sample_id]

# Make a vecor of all common MIMAT before filtering.
ucam_genes = read.table(paste(annotdir, "/", "UCAM_genes.txt", sep=""), header = FALSE, sep="\t", 
												stringsAsFactors=FALSE)[,1] # The IDs of all genes from the UCAM platform, also the ones filtered out in the data matrix.
unfilteredUCAMMIMAT =  mimatmapping[ucam_genes]
rm(ucam_genes)
unfilteredUCAMMIMAT = unfilteredUCAMMIMAT[!is.na(unfilteredUCAMMIMAT)]
filteredAHUSMIMAT =  rownames(ahus_matrix)
unfilteredAHUSMIMAT =  mimatmapping[ahus_uRNAList$genes$GeneName]
unfilteredAHUSMIMAT = unfilteredAHUSMIMAT[!is.na(unfilteredAHUSMIMAT)]
unfilteredcommonMIMAT = intersect(unfilteredAHUSMIMAT, unfilteredUCAMMIMAT)

# # a back to name mapping for the two data sets
# commonanot = data.frame(MIMAT=unfilteredcommonMIMAT, UCAMname=NA, AHUSname=NA)
# commonanot$UCAMname = names(unfilteredUCAMMIMAT[match(unfilteredcommonMIMAT, unfilteredUCAMMIMAT)])
# commonanot$AHUSname = names(unfilteredAHUSMIMAT[match(unfilteredcommonMIMAT, unfilteredAHUSMIMAT)])
# a = (commonanot$UCAMname == commonanot$AHUSname)
# commonanot[!a,] # only 3 differed

# make a mapping back to the microRNA name. only use AHUS names since thay were almost all the same
unfilteredcommonname = names(unfilteredAHUSMIMAT[unfilteredAHUSMIMAT %in% unfilteredcommonMIMAT])
names(unfilteredcommonname) = mimatmapping[unfilteredcommonname]
	


# filter data based on common MIMATS,
filteredcommonMIMAT = intersect(rownames(ucam_matrix), rownames(ahus_matrix))
ucam_matrix = ucam_matrix[filteredcommonMIMAT, ]
ahus_matrix = ahus_matrix[filteredcommonMIMAT, ]
common_matrix= cbind(ahus_matrix, ucam_matrix)

inputisread = TRUE

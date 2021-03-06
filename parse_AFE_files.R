

# setwd("/Users/vegardnygaard/prosjekter/rprojects/Eurocan_microRNA")
# source("parse_AFE_files.R")
starttime = Sys.time()

print( paste("Parsing AFE files started ",as.character(starttime)))
			 
library(AgiMicroRna)
# did not work for all files
# source("/Volumes/abel/prosjekter/eurocan/scripts/common_eurocan.R")

inputdir = "input"
AFE_dir = "not_in_github/AHUS_AFE"
outputdir = "not_in_github"
startdir = getwd()

############################################################
### Read sampleannotation
############################################################
sampleannotation_file="sampleannotation.txt"

sampleannotation=read.table(paste(inputdir,  "/", sampleannotation_file, sep=""),
														sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
														strip.white=TRUE, comment.char ="")
sampleannotation = sampleannotation[sampleannotation$provider=="AHUS",]
row.names(sampleannotation) = sampleannotation$sample_id

############################################################
### Parse AFE files
############################################################
sampleannotation$FileName = sampleannotation$datafile
setwd(AFE_dir) # no way to give a dir to the parser
AHUS_agiMicroRNAdata = readMicroRnaAFE(sampleannotation[,])
setwd(startdir)
############################################################
### Collapse to microRNA with RMA
############################################################
ahus_uRNAList  = rmaMicroRna(AHUS_agiMicroRNAdata, normalize = TRUE, background = FALSE)

############################################################
### Save object for later use
############################################################
save(ahus_uRNAList, file=paste(outputdir, "/AHUS_rma_uRNAList.rdata", sep=""))

#write matrix as text for deposit to arrayExpress or similar
arrayExpressmatrix = ahus_uRNAList$TGS
rownames(arrayExpressmatrix) = ahus_uRNAList$genes$GeneName
write.table(arrayExpressmatrix, file=paste(outputdir, "/toArrayExpress.txt", sep=""),
						row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

		 
sink("parse_AFE_files_sessionInfo.txt")
print( paste("Parsing AFE files ended ",as.character(Sys.time()), "Time spent", 
						 as.integer(round(difftime(Sys.time(),starttime, units="mins"))), "minutes"))
print(sessionInfo())
sink()
		 


		 
		 
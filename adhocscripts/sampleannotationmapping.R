




# simple script to subset a bigger sampleannotation file to the samples and annotation used in this project.
# this is a ad-hoc not runnable notes to self-script


biggerfile = "/Users/vegardnygaard/prosjekter/eurocan/tidy/eurocan_data_summary_2015-06-23/sampleannotation_microRNA.txt"
bigsa=read.table(biggerfile,
														sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
														strip.white=TRUE, comment.char ="")



# the samples used in this project
annotdir = "../input"
sampleannotation_file="sampleannotation.txt"
sa=read.table(paste(annotdir,  "/", sampleannotation_file, sep=""),
														sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
														strip.white=TRUE, comment.char ="")


# attach 
sa$iClust = NA
sa$iClust = bigsa$iClust[ match(sa$fulldatapath, bigsa$fulldatapath) ]
write.table(sa, file=sampleannotation_file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")






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


# attach pam50 scores that where calculated from mRNA samples in a related project.
UCAM_mRNAannto=read.table("../not_in_github/mRNAannot/clinical-data-illumina-with-subtyping.txt",
								 sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
								 strip.white=TRUE, comment.char ="")
All_mRNAannto=read.table("../not_in_github/mRNAannot/sampleannotation_w_subtype_agilent.txt",
													sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE,
													strip.white=TRUE, comment.char ="")

# check that the UCAMs pam50 still are the same
UCAM_mRNAannto$similarid = gsub("MB-", "UCAM_MB.", UCAM_mRNAannto$sample_id)
a = match( UCAM_mRNAannto$similarid, sa$sample_id)
UCAM_mRNAannto = UCAM_mRNAannto[!is.na(a),]
a = match( UCAM_mRNAannto$similarid, sa$sample_id)
UCAM_mRNAannto$subtype.est = gsub("Basal", "Basallike", UCAM_mRNAannto$subtype.est)
UCAM_mRNAannto$subtype.est = gsub("Normal", "Normallike", UCAM_mRNAannto$subtype.est)
table(UCAM_mRNAannto$subtype.est==sa$pam50_est[a]) # ok.
sa$accountedfor=FALSE
sa$accountedfor[a] = TRUE

a =match(All_mRNAannto$sample_id, sa$sample_id)
All_mRNAannto = All_mRNAannto[!is.na(a),]
a =match(All_mRNAannto$sample_id, sa$sample_id)
All_mRNAannto$sub = gsub("^Basal$", "Basallike", All_mRNAannto$sub)
All_mRNAannto$sub = gsub("^Normal$", "Normallike", All_mRNAannto$sub)
b =(!is.na(All_mRNAannto$sub) & All_mRNAannto$sub!=sa$pam50_est[a])
data.frame(All_mRNAannto$sample_id, All_mRNAannto$sub, sa$pam50_est[a], sa$tissue_type[a])[b,]
b =(!is.na(All_mRNAannto$sub) & All_mRNAannto$sub==sa$pam50_est[a])
sa$accountedfor[a][b] = TRUE


sa[!sa$accountedfor, c("sample_id", "pam50_est", "tissue_type", "accountedfor")]


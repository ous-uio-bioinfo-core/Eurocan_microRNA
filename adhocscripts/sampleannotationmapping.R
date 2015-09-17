




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

UCAM_mRNAannto$similarid = gsub("MB-", "UCAM_MB.", UCAM_mRNAannto$sample_id)
a = UCAM_mRNAannto$similarid %in% sa$sample_id
newtab  = data.frame( sample_id=UCAM_mRNAannto$similarid[a], newpam50est=UCAM_mRNAannto$subtype.est[a], stringsAsFactors=FALSE)

a = All_mRNAannto$sample_id %in%  sa$sample_id & All_mRNAannto$provider!="UCAM"
tmptab  = data.frame( sample_id=All_mRNAannto$sample_id[a], newpam50est=All_mRNAannto$sub[a], stringsAsFactors=FALSE)

newtab = rbind(newtab,tmptab)
newtab$newpam50est = gsub("Basal", "Basallike", newtab$newpam50est)
newtab$newpam50est = gsub("Normal", "Normallike", newtab$newpam50est)

comtab = data.frame(sa$sample_id, sa$tissue_type, sa$pam50_est, newpam50=newtab$newpam50est[match(sa$sample_id,newtab$sample_id)], stringsAsFactors=FALSE)
comtab$newpam50[is.na(comtab$newpam50)] = "unknown"
a = comtab$newpam50!=comtab$sa.pam50_est
comtab[a,] # tre var forskjellige.

table(sa$sample_id==comtab$sa.sample_id)
sa$pam50_est = comtab$newpam50
write.table(sa, file=sampleannotation_file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



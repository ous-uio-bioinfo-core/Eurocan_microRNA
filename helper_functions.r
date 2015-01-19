


plotmanypca = function(matrices, annot, plotdirname,  views=c("tissue_provider", "provider", "tissue_type", "platform"))
{  
  
  for(thisname in names(matrices))
  {
    #### be sure to sort sampleannotation as data	
    thismatrix = matrices[[thisname]]
    a = dimnames(thismatrix)[[2]] %in% annot$sample_id
    thismatrix = thismatrix[,a]
    sa = annot[match(dimnames(thismatrix)[[2]] , annot$sample_id), ]
    
    ## handy shortcuts for plot
    sa$tissue_provider = paste(sa$tissue_type, sa$provider, sep="_")
    
    # plot also the pvca plot.
    # plotpvca(matrices[[thisname]], sa, paste(plotdirname, "/pvca_", transcript_type,"_", thisname, ".png", sep=""))
    
    date()
    thisprcomp=prcomp( (t(na.omit(thismatrix))))
    date()
    
    
    for(view in views)
    {
      thiscolorsname=paste(view,"_color", sep="")
      thiscolor=NA
      if(thiscolorsname %in% names(sa))
        thiscolor = sa[,  thiscolorsname]
      plotonepca( name=paste(thisname, view, sep="_"), thisannot=sa[,view], thisplatform=sa[,"platform"] ,thisprcomp=thisprcomp, plotdirname=plotdirname, thiscolors=thiscolor)
    }
  }
}

# name = paste(transcript_type,thisname, view, sep="_")
# thisannot = sa[,view]
# thisplatform = thisplatform=sa[,"platform"]
# thisprcomp 
# plotdirname 

plotonepca = function(name, thisannot, thisplatform, thisprcomp, plotdirname="", thiscolors=NA)
{
  colorpal = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))
  plotsize = 1600
  goodpch = c(19,17,25,23)
  
  if(!is.factor(thisannot))
    thisannot=factor(thisannot)
  
  #print(thisannot)
  if(is.na(thiscolors))
    thiscolors = colorpal[as.numeric(thisannot)]
  
  annotlabels=paste(unique(thisannot),table(thisannot)[unique(thisannot)])
  annotlabelscolors=unique((thiscolors))
  a=match(levels(thisannot), unique(thisannot))
  a=a[!is.na(a)]
  annotlabels=annotlabels[a]
  annotlabelscolors=annotlabelscolors[a]
  
  thispch = factor(thisplatform)
  levels(thispch) =  goodpch

  
  if(plotdirname!="")    
    png(file = paste(plotdirname , "/pca_",name, ".png", sep=""),pointsize = 12, width = plotsize, height = plotsize, units="px")
  orgmar = par()$mar
  orgxpd=par()$xpd
  par(xpd=T, mar=par()$mar+c(0,0,0,6))
  plot(thisprcomp$x,  col=thiscolors, cex=1,  main=name, pch=as.numeric(as.character(thispch)) ,cex.main=1, bg=thiscolors)
  legend("topright" ,inset=c(-0.1,0), legend=annotlabels, text.col="black", cex=1, bty="n", pch=15, col=annotlabelscolors)
  legend("bottomright" ,inset=c(-0.1,0), pch=unique(as.numeric(as.character(thispch))), legend=paste("",unique(thisplatform)), text.col="black", cex=1, bty="n")
  par(mar=orgmar, xpd=orgxpd)
  if(plotdirname!="") 
    dev.off()
  
}




xxxold_makegeneannotation = function(startnames)
{	
  ret = data.frame(Symbol=startnames, stringsAsFactors=FALSE)
  dimnames(ret)[[1]] = startnames
  
  #make the ids comparable to ids from volinia:
  ret$modID1 = gsub("put-" , "" ,  ret$Symbol )
  ret$modID1 = substring( ret$Symbol, regexpr( "-", ret$Symbol)+1 )
  ret$modID2 = gsub("-3p" , "" , ret$modID1 )
  ret$modID2 = gsub("-5p" , "" , ret$modID2 )
  ret$modID2 = gsub("\\*" , "" , ret$modID2 )
  return(ret)
}




# small function to find best probable match between ids.
#thisid="let-7b"
xxxold_getbestprobematch = function(thisid, ga)
{
  a = ga$modID1 == thisid # direct match
  if(sum(a)==1)
  {
    #return(geneannotation$ProbeName[a])
    return(ga$Symbol[a])
  }
  modid = thisid
  modid= gsub("-3p" , "" , modid )
  modid = gsub("-5p" , "" , modid )
  modid = gsub("\\*" , "" , modid )
  modid = gsub("\\(2\\)" , "-2" , modid )
  modid = gsub("-RNASEN" , "" , modid )
  modid = gsub("-DICER1" , "" , modid )
  a = ga$modID2 == modid
  #return(geneannotation$ProbeName[a])
  if(sum(a)==1)
  {
    return(ga$Symbol[a])
  }
  if(sum(a)>1)
  {
    #print(paste("More than one  match for microRNA: ", thisid, "  match:", paste(ga$Symbol[a], collapse=", ")))
    return(ga$Symbol[a])		
  }
  return(vector())
  
}


xxxold_find_diff_genes = function(datamatrix, labels, tissue1, tissue2, mimatnames)
{  

  t1 = labels==tissue1
  t2 = labels==tissue2
  t1[is.na(t1)]=FALSE
  t2[is.na(t2)]=FALSE
  fulltab = data.frame(  dimnames(datamatrix)[[1]],
                         mimatnames[dimnames(datamatrix)[[1]]],
                         rowMeans(datamatrix[,t1], na.rm=TRUE),
                         rowMeans(datamatrix[,t2], na.rm=TRUE))
  names(fulltab) = c("microRNA", "MIMATID", tissue1, tissue2)
  fulltab$totintensity = (fulltab[,tissue1]+ fulltab[,tissue2])/2
  fulltab$fc =  fulltab[,tissue1] - fulltab[,tissue2]  
  fulltab$pval=apply(datamatrix, 1, 
                          function(x) { wilcox.test(x[t1], x[t2])$p.value } )
                     #function(x) { t.test(x[t1], x[t2])$p.value } )
  fulltab$fdr = round(p.adjust(fulltab$pval, "BH"), 3)
  fulltab = fulltab[order(fulltab$pval),]
  return(fulltab)
}

xxxold_compare_with_volinia = function(voliniatab, ourtab, geneannotation, mimatnames, validationFDR=0.05)
{
  
  arraytab= ourtab[, c("logFC", "adj.P.Val")]  
  names(voliniatab)[-2] = paste("Volinia_", names(voliniatab)[-2],  sep="")
  names(arraytab) = paste("array_", names(arraytab),  sep="")
  arraytab$MIMATID = rownames(ourtab)
  comparisontab = merge(voliniatab, arraytab, by.x="MIMATID", by.y="MIMATID", all.x=TRUE)

  comparisontab$directionagreement =  (comparisontab$array_logFC  * comparisontab$Volinia_log2FC) > 0
  comparisontab$confirmed = comparisontab$array_adj.P.Val < validationFDR & comparisontab$directionagreement
  comparisontab = comparisontab[order(comparisontab$array_adj.P.Val), ]
  return(comparisontab)
}

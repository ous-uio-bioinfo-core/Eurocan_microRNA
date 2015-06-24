


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



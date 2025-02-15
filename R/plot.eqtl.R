`plot.eqtl` <- function(x, file=NULL, which=NULL, sig=0.01, verbose=TRUE, centered=TRUE,log=FALSE, mc.cores=1, genome=NULL, exon.annot=NULL, bed.file=NULL, ...){
  if(is.null(x$windowSize)){
    humanGenome68 <- data.frame(chr=c(1:22,"X","Y"),
                                length=c(249250621, 
                                         243199373, 
                                         198022430, 
                                         191154276, 
                                         180915260, 
                                         171115067, 
                                         159138663, 
                                         146364022, 
                                         141213431, 
                                         135534747, 
                                         135006516, 
                                         133851895, 
                                         115169878, 
                                         107349540, 
                                         102531392, 
                                         90354753, 
                                         81195210, 
                                         78077248, 
                                         59128983, 
                                         63025520, 
                                         48129895, 
                                         51304566, 
                                         155270560, 
                                         59373566)) 
    
    if(missing(genome)){
      genome <- humanGenome68
      warning("Warning!!! No genome information provided, use the default (Ensembl Human, build 68).")
    }
    plotTrans(x, genome)
  } else {
    plotSingle.eqtl(x=x, file=file, which=which, sig=sig, verbose=verbose, centered=centered,log=log, exon.annot=exon.annot, bed.file=bed.file, ...)
  }
}

`plotSingle.eqtl` <- function(x, file=NULL, which=NULL, sig=0.01, verbose=TRUE, centered=TRUE,log=FALSE,x2=NULL, annot=NULL, double=FALSE, exon.annot=exon.annot, bed.file=bed.file, ...){
  if(x$type=="full")
  {
    windowSize <- x$windowSize
    x <- x$eqtl
    x2 <- x2$eqtl
    
    if(!is.null(which)==TRUE) x <- x[which]

    if(is.null(file))
    {
      plotIt(x=x,sig=sig,windowSize=windowSize,centered=centered,log=log,x2=x2,annot=annot,double=double, exon.annot=exon.annot, bed.file=bed.file)
    } else {
      pdf(file=file,width=10,height=10)
	plotIt(x=x,sig=sig,windowSize=windowSize,centered=centered,log=log,x2=x2,annot=annot, double=double, exon.annot=exon.annot, bed.file=bed.file)
      dev.off()
    }
    invisible()
  } else {
    if(verbose==TRUE)
    {
      cat("The plot function was applied to the 'sig' eQTL - this means the eQTL has to be performed again for each single plot. \n")
      cat("The 'full' stores all required information to plot it directly, but consumes more memory) \n")
    }
    #function(gex, geno, xAnnot, xSamples=NULL, genoSamples=NULL, windowSize=0.5, method="LM", mc=1, sig=NULL, which=NULL)
    if(!is.null(file))
    {
      pdf(file=file,width=10,height=10)
    }
  
    if(is.null(which)) which <- x$which
  
    for(i in 1:length(which))
    {
      if(sum(is.element(colnames(x$gex),which[i]))>0 & sum(is.element(names(x$xAnnot),which[i]))>0){
	tempEQTL <- eQTL(gex=x$gex,geno=x$geno,xAnnot=x$xAnnot,xSamples=x$xSamples,genoSamples=x$genoSamples,windowSize=x$windowSize,method=x$method,mc=x$mc,which=which[i],sig=NULL, verbose=FALSE)
	windowSize <- tempEQTL$windowSize
	xx <- tempEQTL$eqtl
      	#if(!is.null(which)==TRUE) xx <- xx
	plotIt(x=xx,sig=sig,windowSize=windowSize,centered=centered,log=log,x2=x2,annot=annot, double=double, exon.annot=exon.annot, bed.file=bed.file)
      } else {
	warning("Check your 'which' parameter, there are no expression and/or annotation for the object ",which[i],"\n")
      }
      if(verbose==TRUE)
      {
        cat("We managed to finnish the ",i," plot.\n")
      }
    }
    
    if(!is.null(file))
    {
      dev.off()
    }
    
  }
} 


plotIt <- function(x,sig,windowSize,centered,log,x2,annot, double, exon.annot, bed.file){
  for(gene in 1:length(x))
  {
      temp <- x[[gene]]
      if(!is.null(x2)) temp2 <- x2[[gene]]
      Nlocs <- length(table(temp$ProbeLoc))
      if(Nlocs==0){
	  minX <- as.numeric(temp$GeneInfo[1,][2] - windowSize*10^6)
          maxX <- as.numeric(temp$GeneInfo[1,][3] + windowSize*10^6)

	  if(log==FALSE){
	     if(!is.null(exon.annot)){
	       plotY <- c(-0.2,1.5) 
	     } else {
	       plotY <- c(0,1.5)  
	     }
	     plotYlab <- "p-value"
	  } else {
	     if(!is.null(exon.annot)){
	       plotY <- c(-0.2,5) 
	     } else {
	       plotY <- c(0,5)  
	     }
	     
	     plotYlab <- "-log(p-value)"
	  }

	  if(double==FALSE){
	    plot(c(-10,-10),ylim=plotY,xlim=c(0,20),xlab="Chromosomal Position in MB",ylab=plotYlab,main=paste(names(x)[gene],"- NO SNPs found"),sub=paste("Chr",temp$GeneInfo[1,1],":",prettyNum(minX,big.mark = ","),"-",prettyNum(maxX,big.mark = ",")),yaxt="n",xaxt="n")
	    axis(1,at=seq(0,20,4),labels=seq(round(minX/10^6,1),round(maxX/10^6,1),length.out=6))
	  } else {
	    plot(c(-10,-10),ylim=plotY,xlim=c(0,20),ylab=plotYlab,main=paste(names(x)[gene],"- NO SNPs found"),yaxt="n",xaxt="n")
	  }
	  axis(2,at=c(seq(0,1,0.2),1.1,1.2, 1.3),labels=c(seq(0,1,0.2),"MAF" ,"Mono.","NA"))
	  

	  # Plot the gene position
	  xG <- (mean(unlist(c(temp$GeneInfo[1,][2:3])))-minX)/(windowSize*10^5)
	  lines(c(xG,xG),c(-1,2),lty="dashed")
      
	  # First plot the gene
	  if(!is.null(exon.annot)){
	    tmp.exon <- exon.annot[exon.annot$gene_id==names(x),]
	    tmp.exon$V4 <- (tmp.exon$V4 - minX)/(windowSize*10^5)
	    tmp.exon$V5 <- (tmp.exon$V5 - minX)/(windowSize*10^5)
	    right.side <- 0
	    left.side <- 20
	    for(exonRun in 1:nrow(tmp.exon)){
	      right.side <- max(right.side, tmp.exon$V4, tmp.exon$V5)
	      left.side <- min(left.side, tmp.exon$V4, tmp.exon$V5)
	      rect(tmp.exon$V4, -0.1, tmp.exon$V4, -0.05, col="black")
	    }
	    rect(left.side, -0.1, right.side, -0.05, col=NA, border="blue")
	  }
	  
	  # Now plot the surrounding genes also
	  if(!is.null(exon.annot)){
	    # Extract the correct chromosome
	    chr.exon <- exon.annot[exon.annot$V1==x[[1]]$GeneInfo[["Chr"]],]
	    chr.exon$V4 <- (chr.exon$V4 -minX)/(windowSize*10^5)
	    chr.exon$V5 <- (chr.exon$V5 -minX)/(windowSize*10^5)
	    
	    # Filter for the right location
	    chr.exon <- chr.exon[chr.exon$V4>0,]
	    chr.exon <- chr.exon[chr.exon$V5>0,]
	    chr.exon <- chr.exon[chr.exon$V4<20,]
	    chr.exon <- chr.exon[chr.exon$V5<20,]
	    
	    # Extract the genes of interest
	    local.genes <- unique(chr.exon$gene_id)
	    
	    # And now plot the genes one by one also
	    if(length(local.genes)>0){
	      
	      for(loc_gene_run in 1:length(local.genes)){
	        tmp.chr.exon <- chr.exon[chr.exon$gene_id==local.genes[loc_gene_run],]
	        right.side <- 0
	        left.side <- 20
	        for(exonRun in 1:nrow(tmp.chr.exon)){
	          right.side <- max(right.side, tmp.chr.exon$V4, tmp.chr.exon$V5)
	          left.side <- min(left.side, tmp.chr.exon$V4, tmp.chr.exon$V5)
	          rect(tmp.chr.exon$V4, -0.15, tmp.chr.exon$V4, -0.1, col="black")
	        }
	        rect(left.side, -0.15, right.side, -0.1, col=NA, border="red")
	        
	      }
	    }
	  }
	  
	  
	  } else {
	for(sub in 1:Nlocs)
	{
	  subPos <- temp$ProbeLoc==sub
	  if(centered==TRUE)
	  {
	    minX <- as.numeric(temp$GeneInfo[sub,][2] - windowSize*10^6)
	    maxX <- as.numeric(temp$GeneInfo[sub,][3] + windowSize*10^6)
	  } else {
	    minX <- min(temp$TestedSNP[subPos,4],na.rm=TRUE)
	    maxX <- max(temp$TestedSNP[subPos,4],na.rm=TRUE)
	  }

	  xPos <- (temp$TestedSNP[subPos,4]-minX)/(windowSize*10^5)
	  yPos <- temp$p.values[subPos]
	 
	  bedOut <- data.frame(temp$GeneInfo[1][1,1],
	                       temp$TestedSNP[subPos,4],
	                       temp$TestedSNP[subPos,4],
	                       temp$p.values,
	                       names(x))
	  
	  if(!is.null(bed.file)){
	    write.table(bedOut, file=bed.file, col.names = FALSE, row.names = FALSE, quote=FALSE, sep="\t")
	  }
	  
    if(!is.null(x2)) yPos2 <- temp2$p.values[subPos]
	  # eliminate rounding errors (drop all p-values which are a bit larger than 1 to 1)
	  

          if(log==FALSE){
	     tempPos <- (yPos>1) & (yPos < 1.2)
	     
	     if(!is.null(exon.annot)){
	       plotY <- c(-0.2, 1.5)  
	     } else {
	       plotY <- c(0,1.5)
	     }
	     
	     plotYlab <- "p-value"
	     #yPos[tempPos] <- 1
	     col <- rep("green",length(subPos))
	     col[temp$p.values[subPos] <= sig] <- "red" 
	     col[yPos[subPos] > 1] <- "black"

	  } else {
	     plotYlab <- "-log(p-value)"
	     # Remove the Mono. and NA cases in case of log plots
             tempPos <- (yPos>1)
             yPos[tempPos] <- 1
	     if(!is.null(x2)){
	      tempPos2 <- (yPos2>1)
              yPos2[tempPos2] <- 1
             }
	     # When we use a permutation type of test the result can be zero - add there a very small value 
  	     yPos[yPos==0] <- min(yPos[yPos>0])/10
	     yPos <- -log(yPos,base=10)
	     if(!is.null(x2)){
	        yPos2[yPos2==0] <- min(yPos2[yPos2>0])/10
	        yPos2 <- -log(yPos2,base=10)
	     }
             col <- rep("green",length(subPos))
	     col[yPos > sig] <- "red" 
	     
	     if(!is.null(exon.annot)){
	       plotY <- c(-0.2,max(yPos))  
	     } else {
	       plotY <- c(0,max(yPos))
	     }
	     if(!is.null(x2)){
                col2 <- rep("gold",length(subPos))
   	        col2[yPos2 > sig] <- "steelblue" 
	        plotY2 <- c(0,max(yPos2))
	     }

	  }

	  if(is.null(annot)){
	      if(double==FALSE){
	         plot(c(-10,-10),ylim=plotY,xlim=c(0,20),xlab="Chromosomal Position in MB",ylab=plotYlab,main=paste(names(x)[gene],"-",sub),sub=paste("Chr",temp$GeneInfo[sub,1],":",prettyNum(minX, big.mark = ","),"-",prettyNum(maxX, big.mark = ",") ),yaxt="n",xaxt="n")
                 axis(1,at=seq(0,20,4),labels=seq(round(minX/10^6,1),round(maxX/10^6,1),length.out=6))
              } else {
                 plot(c(-10,-10),ylim=plotY,xlim=c(0,20),ylab=plotYlab,main=paste(names(x)[gene],"-",sub),yaxt="n",xaxt="n")
              }
	  } else {
	      # Put here what happens when we have an annotion (filter them, create the double plot, remove the center line...)
	      ## THIS IS NOT PERFECT YET, I SHOULD USE LAYOUT, SO THAT THE FIELDS CAN BE OF DIFFERENT SIZE!!!
	       
	    oldpar <- par(no.readonly = TRUE)
	    on.exit(par(oldpar))           
	    
	      par(mfrow=c(2,1),
                  oma=c(0,0,2,0),
                  mar=c(0,4.1,0,2.1))
              plot(c(-10,-10),ylim=plotY,xlim=c(0,20),xlab="",ylab=plotYlab,main=paste(names(x)[gene],"-",sub),yaxt="n",xaxt="n")
	  }
	  if(log==FALSE){
	    axis(2,las=2,at=c(seq(0,1,0.2), 1.1, 1.2, 1.3),labels=c(seq(0,1,0.2), "MAF", "Mono.","NA"))
	  } else {
            axis(2,at=c(seq(0,max(yPos),0.5)),labels=c(seq(0,max(yPos),0.5)))
	  }
	  
	  if(log==FALSE){
	     points(xPos,yPos,col=col,pch=20)
	  } else {
	     for(i in 1:length(xPos)){
		if(is.null(x2)){
		   lines(c(xPos[i],xPos[i]),c(0,yPos[i]),col=col[i])
		   if(col[i]=="red") points(xPos[i],yPos[i], col="red",pch=20)
		} else {
		   if(yPos[i]>yPos2[i]){
                       lines(c(xPos[i],xPos[i]),c(0,yPos[i]),col=col[i])
		       if(col[i]=="red") points(xPos[i],yPos[i], col="red",pch=20)
                       lines(c(xPos[i],xPos[i]),c(0,yPos2[i]),col=col2[i])
		       if(col2[i]=="steelblue") points(xPos[i],yPos2[i], col="steelblue",pch=20)

		   } else {
                       lines(c(xPos[i],xPos[i]),c(0,yPos2[i]),col=col2[i])
		       if(col2[i]=="steelblue") points(xPos[i],yPos2[i], col="steelblue",pch=20)
                       lines(c(xPos[i],xPos[i]),c(0,yPos[i]),col=col[i])
		       if(col[i]=="red") points(xPos[i],yPos[i], col="red",pch=20)
		   }
		}
	     }
	  }
	  
	  if(log==FALSE){
	      lines(c(-10,20),c(1,1),lty="dotted")
	  } else {
	      # No other lines at the moment
	  }
	  lines(c(-10,20),c(sig,sig),lty="dotted")

	  # Plot the gene position
	  xG <- (mean(unlist(c(temp$GeneInfo[sub,][2:3])))-minX)/(windowSize*10^5)
	  lines(c(xG,xG),c(-1,2+max(yPos)),lty="dashed")

	  if(!is.null(exon.annot)){
	    tmp.exon <- exon.annot[exon.annot$gene_id==names(x),]
	    tmp.exon$V4 <- (tmp.exon$V4 -minX)/(windowSize*10^5)
	    tmp.exon$V5 <- (tmp.exon$V5 -minX)/(windowSize*10^5)
	    right.side <- 0
	    left.side <- 20
	    for(exonRun in 1:nrow(tmp.exon)){
	      right.side <- max(right.side, tmp.exon$V4, tmp.exon$V5)
	      left.side <- min(left.side, tmp.exon$V4, tmp.exon$V5)
	      rect(tmp.exon$V4, -0.1, tmp.exon$V4, -0.05, col="black")
	    }
	    rect(left.side, -0.1, right.side, -0.05, col=NA, border="blue")
	  }
	  
	  # Then plot the surrounding genes
	  if(!is.null(exon.annot)){
	    # Extract the correct chromosome
	    chr.exon <- exon.annot[exon.annot$V1==x[[1]]$GeneInfo[["Chr"]],]
	    chr.exon$V4 <- (chr.exon$V4 -minX)/(windowSize*10^5)
	    chr.exon$V5 <- (chr.exon$V5 -minX)/(windowSize*10^5)
	   
	    # Filter for the right location
	    chr.exon <- chr.exon[chr.exon$V4>0,]
	    chr.exon <- chr.exon[chr.exon$V5>0,]
	    chr.exon <- chr.exon[chr.exon$V4<20,]
	    chr.exon <- chr.exon[chr.exon$V5<20,]
	     
	    # Extract the genes of interest
	    local.genes <- unique(chr.exon$gene_id)
	    
	    # And now plot the genes one by one also
	    if(length(local.genes)>0){
	      
	      for(loc_gene_run in 1:length(local.genes)){
	        tmp.chr.exon <- chr.exon[chr.exon$gene_id==local.genes[loc_gene_run],]
	        right.side <- 0
	        left.side <- 20
	        for(exonRun in 1:nrow(tmp.chr.exon)){
	          right.side <- max(right.side, tmp.chr.exon$V4, tmp.chr.exon$V5)
	          left.side <- min(left.side, tmp.chr.exon$V4, tmp.chr.exon$V5)
	          rect(tmp.chr.exon$V4, -0.15, tmp.chr.exon$V4, -0.1, col="black")
	        }
	        rect(left.side, -0.15, right.side, -0.1, col=NA, border="red")
	        
	      }
	    }
	  }
	  
	  
	  # Now plot the annotation spot
	  if(!is.null(annot)){
     	        oldpar <- par(no.readonly = TRUE)
	            on.exit(par(oldpar))         
	    
              par(oma=c(1,0,0,0),
                  mar=c(3,4.1,0.1,2.1), new=TRUE)
	      plot(c(-10,-10),ylim=c(0,1),xlim=c(0,20),xlab="Chromosomal Position in MB",ylab="",main="",sub=paste("Chr",temp$GeneInfo[sub,1],":",minX,"-",maxX),xaxt="n")
              axis(1,at=seq(0,20,4),labels=seq(round(minX/10^6,1),round(maxX/10^6,1),length.out=6))
	      for(annotRun in 1:nrow(annot)){
		lines(c(annot[annotRun,3]/10^6,annot[annotRun,4]/10^6),c(0.5,0.5),lwd=2)
		lines(c(40,60),c(0.5,0.5))
	      }
	  }
	}
      }
  }
}

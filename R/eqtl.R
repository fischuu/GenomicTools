eQTL <- function(gex, geno, xAnnot=NULL, xSamples=NULL, genoSamples=NULL, windowSize=0.5, method="LM", mc=1, sig=NULL, which=NULL,nper=2000 , usehoardeR=FALSE, verbose=TRUE){

    if(usehoardeR) stop("Currently not support, please keep the default setting!")
  
  # If the annotations are given as data frame, we will transform them into a list
    if(is.data.frame(xAnnot)){
      if(is.factor(xAnnot[,1])) xAnnot[,1] <- as.character(xAnnot[,1])
      if(is.factor(xAnnot[,2])) xAnnot[,2] <- as.character(xAnnot[,2])
      if(verbose==TRUE) cat("We will transform the gene annotations into a list (",date(),")!\n", sep="")
      xAnnot <- makeAnnotList(xAnnot)
      if(verbose==TRUE) cat("We transformed the gene annotations (",date(),")!\n", sep="")
    }

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
    if(is.character(geno)==TRUE)
    {
      # Check if there is a ped/map pair or vcf given as input
      # CASE: VCF
        if(toupper(substr(geno,nchar(geno)-2, nchar(geno)))=="VCF"){
          stop("The vcf format is not yet supported. Please use the ped/map format.")
      # CASE: PED/MAP
        } else {
          if(verbose==TRUE) cat("Start reading the genotype information at",date(),"(ped/map format assumed)\n")
            if(usehoardeR){
              genotData <- importPedMap(ped=geno)            
            } else {
              genotData <- read.pedfile(file=paste(geno,".ped",sep=""), snps=paste(geno,".map",sep=""))
            }
        }
    # Case: No string provided, assume that genotype data was read in properly
      } else {
      # XXX: Include here still a input check!!!
        genotData <- geno
      }

  # Input checks
    if(is.vector(gex) & is.null(xAnnot))
    {
      warning("No annotations given, we will test all given SNPs against all given expressions!\n")
      xAnnot <- data.frame(Gene=as.character(genotData$map[,2]),Chr=as.character(genotData$map[,1]),Start=genotData$map[,4],End=genotData$map[,4])
      xAnnot <- makeAnnotList(xAnnot)
      windowSize=NULL
    }
    if(is.matrix(gex) & is.null(xAnnot))
    {
      warning("No annotations given, we will test all given SNPs against all given expressions!\n gex is a matrix, hence this might take a while...\n")
      xAnnot <- data.frame(gene="",Chr=as.character(genotData$map[1,1]),Start=genotData$map[1,4],End=genotData$map[1,4])
  #    xAnnot <- makeAnnotList(xAnnot)
      windowSize=NULL
      
    }

  # If no separate genoSamples object is given, we take those from the snpStats object:
    if(is.null(genoSamples)) genoSamples <-  as.character(genotData$fam$member)
  
  # Take only those genes from the annotation list that were requested
    if(!is.null(which)) {
      xAnnot <- xAnnot[is.element(names(xAnnot),which)]
      gex <- gex[is.element(names(gex),which)]
    }
  
  # Input checks
    single <- FALSE
    th <- windowSize
    method <- match.arg(method,c("LM","directional"))

  # If only one gene is given it could be a vector, transform it here to a column matrix
    if(is.vector(gex)){
      tempNames <- names(gex)
      if(length(xAnnot)>1){
        warning("Vector of gene expression is provided together with several annotations. Gene expressions will be repeated for every annotation!")
        cat("Vector of gene expression is provided together with several annotations. Gene expressions will be repeated for every annotation!\n")
        gex <- matrix(rep(temp,length(xAnnot)),ncol=length(xAnnot), byrow=FALSE)
      } else {
        gex <- matrix(gex,ncol=1)
        single <- TRUE
      }
        if(!is.null(xSamples)) tempNames <- xSamples
        if(is.null(tempNames)){
          cat("Vector of expression values in unnamed. We assume same order in expression and genotype objects and match samples based on that.\n")
          tempNames <- as.character((genoData$fam)[,1])
        }
        rownames(gex) <- tempNames
        colnames(gex) <- names(xAnnot)
        gexColNames <- colnames(gex)
        
    } else {
      if(!is.null(xSamples)){
        if(nrow(gex)!=length(xSamples)){
          stop("nrow(gex)!=length(xSamples): The number of names given in xSample has to be the same as needed for the expression matrix!")
        } else {
          rownames(gex) <- xSamples  
        }
      } 
    }
 
  # In case that the row names have been changed, bring them into an order
    rownames(genotData$map) <- 1:nrow(genotData$map)

  # Sample statistics
    overlap <- is.element(rownames(gex),genoSamples)
    olPerc <- round(sum(overlap)/nrow(gex)*100,1)
    if(sum(overlap)==0) stop("No matching expression / genotype sample names!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the samples in the expression data the genotype information available. \n")

  # Location statistics
    overlap <- is.element(colnames(gex),names(xAnnot))
    olPerc <- round(sum(overlap)/ncol(gex)*100,1)
    if(sum(overlap)==0) stop("No matching expression probe names / probe name annotations!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the expression data the annotations. \n")

  # Probe statistics
    matchingProbes <- colnames(gex)[is.element(colnames(gex),names(xAnnot))]
    if(verbose==TRUE) cat("We will investigate for",length(matchingProbes),"genes possible eQTLs! \n")
    result <- list()

  # Reducing the expression data to those rows, where we have also genotype information available
    gex <- gex[is.element(rownames(gex),genoSamples),]
    if(single==TRUE){
      gex <- t(t(gex))
      colnames(gex) <- gexColNames
    }

  # Now go through all possible probes
    eqtl <- list()
 
    for(probeRun in 1:length(matchingProbes)){
    # Do that for each possible location of the probe (might not be unique...)
      tempAnnot <- xAnnot[[which((names(xAnnot)==matchingProbes[probeRun])==TRUE)]]
      eqtlTemp <- list()
   
      for(tempRun in 1:nrow(tempAnnot)){
      # SNP locations of variants inside the provided window
        SNPloc <- getSNPlocations(genotInfo=genotData$map,annot=tempAnnot[tempRun,],th=th)
  
      # Run eQTL only if there are SNPs inside the window
        if(dim(SNPloc$SNPloc)[1]>0){
          if(usehoardeR){
            SNPmatrix <- genotData$geno[,SNPloc$SNPcol]
            genoGroups <- rearrange(SNPmatrix,rownames(gex),genoSamples)
          } else {
            SNPmatrix <- genotData$genotypes[,SNPloc$SNPcol]
            genoGroups <- as(SNPmatrix,"numeric")
            genoGroups <- rearrange(genoGroups,rownames(gex),genoSamples)
          }

        # eQTL case : LM
          if(method=="LM"){
          # if sig is set to Null all results will be reported - This might be very memory consuming!!!
  	        if(is.null(sig)){
  	            eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlLM(genoGroups,gex[,probeRun], mc=mc))
  	        } else {
  	            p.values <- eqtlLM(genoGroups,gex[,probeRun], mc=mc)
            	  pPos <- p.values<=sig
  	            eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
  	        }
        # eQTL case: directional
          } else if(method=="directional"){
          # if sig is set to Null all results will be reported - This might be very memory consuming!!!
            if(is.null(sig)){ 
  	           eqtlTemp[[tempRun]] <- list(ProbeLoc=rep(tempRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlDir.P(genoGroups,gex[,probeRun],mc=mc,nper=nper))
  	        } else {
  	           p.values <- eqtlDir.P(genoGroups,gex[,probeRun],mc=mc,nper=nper)
  	           pPos <- p.values<=sig
  	           eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
          	}
          } 
       } else {
          warning("There were no variants within the window of gene ", matchingProbes[probeRun])
          p.values <- 2  # Take a p-value of 2 as internal indicator that there wasn't anything to test!
          eqtlTemp[[tempRun]] <- data.frame(ProbeLoc=-1,TestedSNP=-1,p.values=-1, Assoc.Gene="-1")
       }
      }
      
      # Join the output
      if(is.null(sig))
      {
        eqtl[[probeRun]] <- joinEQTL(eqtlTemp)
        eqtl[[probeRun]]$GeneInfo <- tempAnnot
      } else {
        if(sum(p.values<=sig)==0){
          eqtl[[probeRun]] <-  data.frame(chr="-1", SNP="-1", Location="-1", p.value="-1", Assoc.Gene=matchingProbes[probeRun], stringsAsFactors=FALSE)
        } else {
          #bedTemp <- joinEQTLsig(eqtlTemp)
          bedTemp <- do.call("rbind", eqtlTemp)
          bedTemp <- cbind(bedTemp,rep(matchingProbes[probeRun],max(1,nrow(bedTemp))))
          colnames(bedTemp) <- c("chr", "SNP", "Location", "p.value", "Assoc.Gene")
          eqtl[[probeRun]] <- bedTemp          
        }
      }
      if(verbose==TRUE) cat ("We calculated eQTLs for gene ",matchingProbes[probeRun]," for ",length(eqtl[[probeRun]]$p.values)," SNPs (",date(),")\n", sep="")
    }

  # Return the result
  if(is.null(sig))
  {
    names(eqtl) <- matchingProbes
    result <- list(bed=NULL,eqtl=eqtl,gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="full")
  } else {
    resBed <- do.call("rbind", eqtl)
    resBed <- resBed[-which(resBed[,1]== -1),]
    result <- list(bed= resBed,gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="sig")
  }
  class(result) <- "eqtl"
  result
}




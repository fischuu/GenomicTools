`print.qtlRes` <- function(x, which=NULL, sig=0.001, ...){
   
  if(is.null(which)) which <- 1:length(x$qtl)
  
  for(phenoRun in which){
    tempQTL <- x$qtl[phenoRun][[1]]
  
    takeThese <- which(tempQTL$p.values<=sig)
    locNames <- NA
    if(length(takeThese)>0){
      tmpLocs <- tempQTL$TestedSNP[takeThese,]
      locNames <<- colnames(tmpLocs)
      tmpP <- tempQTL$p.values[takeThese]
      phenoOut <- cbind(tmpLocs,tmpP, colnames(x$pheno)[phenoRun])
      ifelse(exists("output"), output <- rbind(output,phenoOut), output <- phenoOut)
    }
  }
  if(!exists("output")){
    output <- c(NA,NA, NA, NA, NA)
    locname <- c("Chr", "Pos")
    warning("No significant results with p <",sig)
  }
  colnames(output) <- c(locNames, "p.value", "assoc. pheno.")
  print(output)
} 



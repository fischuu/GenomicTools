`print.qtlRes` <- function(x, which=NULL, sig=0.001, ...){
   
  if(is.null(which)) which <- 1:length(x$qtl)
  
  for(phenoRun in which){
    tempQTL <- x$qtl[phenoRun][[1]]
  
    takeThese <- which(tempQTL$p.values<=sig)
      
    tmpLocs <- tempQTL$TestedSNP[takeThese,]
    tmpP <- tempQTL$p.values[takeThese]
    phenoOut <- cbind(tmpLocs,tmpP, colnames(x$pheno)[phenoRun])
    if(phenoRun==1) output <- phenoOut
    if(phenoRun>1) output <- rbind(output,phenoOut)
  }
  colnames(output) <- c("chr","snp.name"," ", "location", "allele.1", "allele.2", "p.value", "assoc. pheno.")
  output
} 


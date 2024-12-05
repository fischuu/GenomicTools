eqtlDir <- function(genoGroups, gex, mc=1, nper, testType, MAF=0){

    innerFunction <- function(k)
    {
      # Calculate MAF
      maf.snp <- sum(genoGroups[genoGroups[,k]<3,k], na.rm=TRUE)/(2*(sum(genoGroups[,k]<3, na.rm=TRUE) ))
      maf.snp <- min(1-maf.snp, maf.snp)
      if(is.nan(maf.snp)) maf.snp <- 0
      
      # Genotypes are coded as '3' if there is missing data, only proced if there is valid data
      if(sum(genoGroups[,k]==3)<(nrow(genoGroups))){
        missingData <- which((genoGroups[,k]==3)==TRUE)
      	PI <- length(table(genoGroups[genoGroups[,k]!=3,k]))
	    
      # All individuals have the same genotype, do nothing
	      if(PI==1){
	        innerOut <- 1.2
	
	      # Then the 2 groups comparison
	      } else if (PI==2){
	        if(testType=="asymptotic") testType <- "external"
	        innerOut <- gmw(gex[genoGroups[,k]!=3],genoGroups[genoGroups[,k]!=3,k],test="mw",type=testType,alternative="two.sided",nper=nper)$p.values

	      # And the three group comparison
	      } else if (PI==3){
	         p1 <- gmw(gex[genoGroups[,k]!=3],genoGroups[genoGroups[,k]!=3,k],test="triple",type=testType,alternative="greater",nper=nper)$p.values
	         p2 <- gmw(gex[genoGroups[,k]!=3],createGroups(genoGroups[genoGroups[,k]!=3,k],c(2,1,0)),test="triple",type=testType,alternative="greater",nper=nper)$p.values
	         innerOut <- min(2*min(p1,p2),1)
	  
	      }
      	
      	# Compute medians for each group
      	groupMedians <- tapply(gex[genoGroups[,k]!=3], genoGroups[genoGroups[,k]!=3,k], median, na.rm=TRUE)
      	
      	
      } else {
        innerOut <- 1.3
        groupMedians <- NA
      }
      
      if(maf.snp<=MAF){
        innerOut <- 1.1
        groupMedians <- NA
      }   
      list(p_value = innerOut, medians = groupMedians)
    }
    
    if (is.matrix(genoGroups)) {
      output <- mclapply(1:ncol(genoGroups), innerFunction, mc.cores=mc)
    } else {
      genoGroups <- as.matrix(genoGroups)
      if (ncol(genoGroups) > 1) genoGroups <- t(genoGroups)
      output <- list(innerFunction(1))
    }
    
    output
}

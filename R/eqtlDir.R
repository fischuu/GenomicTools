eqtlDir.P <- function(genoGroups,gex,mc=mc,nper){

  output <- c()
    innerFunction <- function(i)
    {
      # Genotypes are coded as '3' if there is missing data, only proced if there is valid data
      if(sum(genoGroups[,i]==3)<(nrow(genoGroups))){
        missingData <- which((genoGroups[,i]==3)==TRUE)
	PI <- length(table(genoGroups[genoGroups[,i]!=3,i]))
	# All individuals have the same genotype, do nothing
	if(PI==1){
	  output[i] <- 1.25
	
	# Then the 2 groups comparison
	} else if (PI==2){
	  output[i] <- gmw(gex[genoGroups[,i]!=3],genoGroups[genoGroups[,i]!=3,i],test="mw",type="permut",alternative="two.sided",nper=nper)$p.values

	# And the three group comparison
	} else if (PI==3){
	  p1 <- gmw(gex[genoGroups[,i]!=3],genoGroups[genoGroups[,i]!=3,i],test="triple",type="permutation",alternative="greater",nper=nper)$p.values
	  p2 <- gmw(gex[genoGroups[,i]!=3],createGroups(genoGroups[genoGroups[,i]!=3,i],c(2,1,0)),test="triple",type="permutation",alternative="greater",nper=nper)$p.values
	  output[i] <- min(2*min(p1,p2),1)
	  
	}
      } else {
	output[i] <- 1.5
      }
   }
  output <- unlist(mclapply(1:ncol(genoGroups),innerFunction,mc.cores=mc))
  output
}

eqtlLM.internal <- function(geno,gex){
    res <- 1.5
    
    if(sum(is.na(geno)) < length(geno)){
      temp <- lm(gex~geno)
      ifelse(is.na(temp$coefficients[2]), res <- 1.25 ,  res <- summary(temp)$coefficients[2,4])
    } 
}

eqtlLM <- function(genoGroups, gex, mc=mc){
  if(mc==1){
    res <- as.vector(apply(genoGroups,2,eqtlLM.internal, gex=gex))
  } else {
    res <- as.vector(unlist(mclapply(1:dim(genoGroups)[2], function(i) eqtlLM.internal(geno=genoGroups[,i], gex=gex), mc.cores=mc)))
  }
  res
}

eqtlLM.naive <- function(genoGroups,gex){

  output <- c()
  for(i in 1:ncol(genoGroups))
  {
    if(sum(is.na(genoGroups[,i]))<(nrow(genoGroups))){
      fitData <- data.frame(x=genoGroups[,i],y=gex)
      temp <- lm(y~x,data=fitData)
      ifelse(is.na(temp$coefficients[2]), output[i] <- 1.25 ,  output[i] <- summary(temp)$coefficients[2,4])
    } else {
      output[i] <- 1.5
    }
  }
  output
}
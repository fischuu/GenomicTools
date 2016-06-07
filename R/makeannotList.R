makeAnnotList <- function(xAnnot){
  result <- list()
  if(sum(colnames(xAnnot)==c("Gene","Chr","Start","End"))!=4){
   stop("Please relabel the column names of xAnnot: 'Gene','Chr','Start','End'")
  }
  for(i in 1:nrow(xAnnot)){
    result[[i]] <- xAnnot[i,2:4] 
  }
  names(result) <- xAnnot[,1]
  result
}
summary.mdr <- function(object, ...){
  
  mdr <- object$mdr

  cat("MDR Summary\n")
  cat("---------------\n")
  cat("Order of SNP models     :",object$fold,"\n")
  cat("# of tested variables   :",dim(object$X)[2],"\n")
  cat("# of used observations  :",dim(object$X)[1],"\n")
  
  invisible(object)
} 

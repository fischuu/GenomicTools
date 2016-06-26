summary.mdr <- function(object, ...){
  
  mdr <- object$mdr

  cat("MDR Summary\n")
  cat("---------------\n")
  cat("Order of SNP models  :",object$fold,"\n")
  invisible(object)
} 

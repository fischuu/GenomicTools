summary.PedMap <- function(object, ...){
  
  cat("PedMap Summary\n")
  cat("---------------\n")
  cat("# of imported SNPs       :",nrow(object$map),"\n")
  cat("# of imported samples    :",nrow(object$fam),"\n")
  cat("# of missing sites       :",object$meta$missing,"\n")
  cat("# of monomorphic sites   :",object$meta$mono,"\n")
  cat("# of multiallelic sites  :",object$meta$multiallelic,"\n")
  cat("Used ped file            :",object$meta$pedFile,"\n")
  invisible(object)
} 

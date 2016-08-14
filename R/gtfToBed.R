gtfToBed <- function(gtf){
  if(sum(class(gtf)=="gtf")==0) stop("This function requires a gtf-object as import (as given by the importGTF() function)")
 bedOut <- data.frame(Chr=gtf$V1,
                      Start=gtf$V4,
                      End=gtf$V5,
                      Gene=gtf$gene_id,
                      stringsAsFactors=FALSE)
 class(bedOut) <- append(class(bedOut), "bed")
 bedOut
}
gtfToBed <- function(gtf){
  if(sum(class(gtf)=="gtf")==0) stop("This function requires a gtf-object as import (as given by the importGTF() function)")
 bedOut <- data.frame(Chr=gtf[,1],
                      Start=gtf[,4],
                      End=gtf[,5],
                      Gene=gtf$gene_id,
                      stringsAsFactors=FALSE)

 colnames(bedOut) <- c("Chr", "Start", "End", "Gene")
 class(bedOut) <- append(class(bedOut), "bed")
 bedOut
}
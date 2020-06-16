x <- "ATTGCGAGCT"

revcomp <- function(x, rev=TRUE){
  out <- unlist(strsplit(x,""))
  for(i in 1:length(out)){
    if(out[i]=="A"){
      out[i] <- "T"
    }  else if(out[i]=="T"){
       out[i] <- "A"
      } else if(out[i]=="C"){
        out[i] <- "G"
      } else if(out[i]=="G") {
        out[i] <- "C"
        }
  }
 
  if(rev){
    out <- paste(out[(length(out):1)], collapse="") 
  } else {
    out <- paste(out, collapse="")
  }
  
  out
  
}
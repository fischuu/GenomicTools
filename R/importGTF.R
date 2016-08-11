importGTF <- function(gtf, skip=0, nrow=-1, use.data.table=FALSE){
  if(use.data.table){
    out <- fread(input=gtf, sep="\t", header=FALSE, skip=skip, nrows=nrow)
  } else {
    
    # Load the file and split the Variable V9
      cuffLoaded <- read.csv(file=gtf, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=skip, nrow=nrow)
      V9 <- cuffLoaded$V9
      V9 <- strsplit(V9,"; ")
      
    # Now get the part before the space (=first) and after (=second)
      first <- sapply(V9, function(x)unlist(regmatches( x , gregexpr( "[_a-zA-Z]+ " , x ) )))
      second <- sapply(V9, function(x)unlist(regmatches( x , gregexpr( "[_a-zA-Z0-9=.]+" , x ) )))
        
    # Column names
      tags <- unique( unlist(first) )
      
    # Intermediate matrices
      temp <- mapply( cbind , second , first )
      
    # Match to appropriate columns and coerce to data.frame
      out <- data.frame( do.call( rbind , lapply( temp , function(x) x[ match( tags , x[,2] ) ]  ) ) , stringsAsFactors=FALSE)
      names(out) <- tags
      out <- cbind(cuffLoaded[,-9],out)
  }
  out
}
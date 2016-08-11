importGTF <- function(file, skip=0, nrow=-1, use.data.table=FALSE, level="gene"){
  gtf <- file
  if(use.data.table){
    cuffLoaded <- fread(input = paste('zcat',gtf), skip=5, colClasses = c("character",
                                                                          "character",
                                                                          "character",
                                                                          "integer",
                                                                          "integer",
                                                                          "character",
                                                                          "character",
                                                                          "character",
                                                                          "character"))
    if(!is.null(level)) cuffLoaded <- cuffLoaded[cuffLoaded$V3==level,]
    if(nrow>0) cuffLoaded <- cuffLoaded[1:nrow,]
  } else {
    cuffLoaded <- read.csv(file=gtf, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=skip, nrow=nrow)
  }
  # Split the variable V9
    V9 <- cuffLoaded$V9
    V9 <- strsplit(V9,"; ")
      
  # Now split the data
    V9Splitted <- sapply(V9, function(x)unlist(x))
        
  # Column names
    tags <- unique( unlist(first) )
      
  # Intermediate matrices
    temp <- mapply( cbind , second , first )
      
  # Match to appropriate columns and coerce to data.frame
    out <- data.frame( do.call( rbind , lapply( temp , function(x) x[ match( tags , x[,2] ) ]  ) ) , stringsAsFactors=FALSE)
    names(out) <- tags
    out <- cbind(cuffLoaded[,-9],out)

    out
}
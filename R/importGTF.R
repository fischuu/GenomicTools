importGTF <- function(file, skip=0, nrow=-1, use.data.table=TRUE, level="gene"){
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
    stop("Currently the importGTF function supports only data tables.")
    cuffLoaded <- read.csv(file=gtf, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=skip, nrow=nrow)
  }
  # Split the variable V9
    V9 <- cuffLoaded$V9
    V9 <- strsplit(V9,"; ")
      
  # Now get the required information from V9
    V9.gene_id <- sapply(V9, function(x) x[grepl("gene_id",x)])
    V9.gene_name <- sapply(V9, function(x) x[grepl("gene_name",x)])
    V9.gene_biotype <- sapply(V9, function(x) x[grepl("gene_biotype",x)])
 
  # Remove the non-informative aprts from that vectors
    V9.gene_id <- gsub("gene_id ","",V9.gene_id)
    V9.gene_name <- gsub("gene_name ","",V9.gene_name)
    V9.gene_biotype <- gsub("gene_biotype ","",V9.gene_biotype)
    V9.gene_id <- gsub('\"',"",V9.gene_id)
    V9.gene_name <- gsub('\"',"",V9.gene_name)
    V9.gene_biotype <- gsub('\"',"",V9.gene_biotype)
    
    cuffLoaded[,V9:=NULL]
    cuffLoaded[,gene_id:=V9.gene_id]
    cuffLoaded[,gene_name:=V9.gene_name]
    cuffLoaded[,gene_biotype:=V9.gene_biotype]
     
    class(cuffLoaded) <- append(class(cuffLoaded), "gtf")
    cuffLoaded
}
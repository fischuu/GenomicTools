importGTF <- function(file, skip=5, nrow=-1, use.data.table=TRUE, level="gene", features=NULL, print.features=FALSE){
  gtf <- file
  if(use.data.table){
    if(last(strsplit(gtf,"\\.")[[1]])=="gz"){
      cuffLoaded <- fread(input = paste('zcat',gtf), skip=skip, colClasses = c("character",
                                                                               "character",
                                                                               "character",
                                                                               "integer",
                                                                               "integer",
                                                                               "character",
                                                                               "character",
                                                                               "character",
                                                                               "character"))
    } else {
      cuffLoaded <- fread(input = gtf, skip=skip, colClasses = c("character",
                                                                 "character",
                                                                 "character",
                                                                 "integer",
                                                                 "integer",
                                                                 "character",
                                                                 "character",
                                                                 "character",
                                                                 "character"))      
    }

    if(!is.null(level)) cuffLoaded <- cuffLoaded[cuffLoaded$V3==level,]
    if(nrow>0) cuffLoaded <- cuffLoaded[1:nrow,]
  } else {
    stop("Currently the importGTF function supports only data tables.")
    cuffLoaded <- read.csv(file=gtf, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=skip, nrow=nrow)
  }
  # Split the variable V9
    V9 <- cuffLoaded$V9
    V9 <- strsplit(V9,"; ")

  # Print the features, if requested
    if(print.features){
      cat("List of features in column 9:\n")
      cat("-----------------------------\n")
      cat(gsub(" ","",paste(sapply(strsplit(V9[[1]]," "),"[",1),"\n")))
    }
    
  # Remove the non-informative aprts from that vectors
    if(is.null(features)){
    # Now get the required information from V9
      gene_id <- sapply(V9, function(x) x[grepl("gene_id",x)])
      gene_name <- sapply(V9, function(x) x[grepl("gene_name",x)])
      gene_biotype <- sapply(V9, function(x) x[grepl("gene_biotype",x)])
      gene_id <- gsub("gene_id ","",gene_id)
      gene_name <- gsub("gene_name ","",gene_name)
      gene_biotype <- gsub("gene_biotype ","",gene_biotype)
      gene_id <- gsub('\"',"",gene_id)
      gene_name <- gsub('\"',"",gene_name)
      gene_biotype <- gsub('\"',"",gene_biotype)
      
      cuffLoaded[,V9:=NULL]
      cuffLoaded[,gene_id:=gene_id]
      cuffLoaded[,gene_name:=gene_name]
      cuffLoaded[,gene_biotype:=gene_biotype]
      
    } else {
      for(frun in 1:length(features)){
        tmpFeature <- sapply(V9, function(x) x[grepl(features[frun],x)])
        tmpFeature <- gsub(" ","",tmpFeature)
        tmpFeature <- gsub(eval(features[frun]),"",tmpFeature)
        tmpFeature <- gsub('\"',"",tmpFeature)
        cuffLoaded[,eval(features[frun]):=tmpFeature]
      }
      cuffLoaded[,V9:=NULL]
    }

    

     
    class(cuffLoaded) <- append(class(cuffLoaded), "gtf")
    cuffLoaded
}
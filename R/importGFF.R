importGFF.internal <- function(file, skip=auto, nrow=-1, use.data.table=TRUE, level="gene", features=NULL, num.features=num.features, print.features=FALSE, merge.feature=NULL, verbose=FALSE){
  gff <- file

   if(skip=="auto"){
     if(last(strsplit(gff,"\\.")[[1]])=="gz"){
       cat("dfg")
     } else {
       con  <- file(gff, open = "r")
       search <- TRUE
       obsRow <- 0
       while(search){
         obsRow <- obsRow + 1
         tmp <- readLines(con, n = 1, warn = FALSE)  
         if(substr(tmp,1,1)!="#"){
           skip <- obsRow -1
           search <- FALSE
         }
         if(obsRow==1000) search <- FALSE
       }
    close(con)
     }
   } else {
     if(!is.numeric(skip)) stop("ERROR: skip needs to be a numeric value!")
   }

  if(verbose) cat("Automatically detected number of rows to skip: ", skip,"\n")
  
  if(use.data.table){
    if(last(strsplit(gff,"\\.")[[1]])=="gz"){
      cuffLoaded <- fread(input = paste('zcat',gff), sep ="\t", skip=skip, colClasses = c("character",
                                                                                          "character",
                                                                                          "character",
                                                                                          "integer",
                                                                                          "integer",
                                                                                          "character",
                                                                                          "character",
                                                                                          "character",
                                                                                          "character"))
              } else {
      cuffLoaded <- fread(input = gff, sep="\t", skip=skip, colClasses = c("character",
                                                                           "character",
                                                                           "character",
                                                                           "integer",
                                                                           "integer",
                                                                           "character",
                                                                           "character",
                                                                           "character",
                                                                           "character"))      
    }

    if(verbose){
      if(sum(cuffLoaded$V3==level)==0){
        availLevels <- unique(cuffLoaded$V3)
        stop("The given import level '",level,"' was not found in the gff! Available level are: \n", paste(availLevels, collapse=" \n "))
      }
    }
    
    if(!is.null(level)) cuffLoaded <- cuffLoaded[cuffLoaded$V3==level,]
    if(nrow>0) cuffLoaded <- cuffLoaded[1:nrow,]
  } else {
    stop("Currently the importGFF function supports only data tables.")
    cuffLoaded <- read.csv(file=gff, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=skip, nrow=nrow)
  }
  # Split the variable V9
    V9 <- cuffLoaded$V9
    V9 <- strsplit(V9,";")

  # Print the features, if requested
    if(print.features || verbose){
      cat("List of features in column 9:\n")
      cat("-----------------------------\n")
      cat(paste(paste(sapply(strsplit(V9[[1]],"="),"[",1),"\n"), collapse=""))
    }
    
  # Remove the non-informative aprts from that vectors
    if(is.null(features)){
    # Now get the required information from V9
      gene_id <- sapply(sapply(V9, function(x) x[grepl("ID",x)]),"[",1)
      gene_name <- sapply(V9, function(x) x[grepl("Name",x)])
      gene_biotype <- sapply(V9, function(x) x[grepl("gene_biotype",x)])
      gene_id <- gsub("ID=","",gene_id)
      gene_name <- gsub("Name=","",gene_name)
      gene_biotype <- gsub("gene_biotype=","",gene_biotype)
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
        tmpFeature <- gsub(";","",tmpFeature)
        tmpFeature <- gsub(eval(features[frun]),"",tmpFeature)
        tmpFeature <- gsub('\"',"",tmpFeature)
        if(sum(is.element(features[frun],num.features))>0) tmpFeature <- as.numeric(tmpFeature)
        cuffLoaded[,eval(features[frun]):=tmpFeature]
      }
      cuffLoaded[,V9:=NULL]
    }

    class(cuffLoaded) <- append(class(cuffLoaded), "gff")
    cuffLoaded
}

importGFF <- function(file, skip="auto", nrow=-1, use.data.table=TRUE, level="gene", features=NULL, num.features=c("FPKM", "TPM"), print.features=FALSE, merge.feature=NULL, merge.all=TRUE, class.names=NULL, verbose=TRUE){

# If no merge feature is given, we assume that only a single gff is to be imported
   if(is.null(merge.feature)){
   out <- importGFF.internal(file=file, skip=skip, nrow=nrow, use.data.table=use.data.table, level=level, features=features, num.features=num.features, print.features=print.features, verbose=verbose)
 } else {
# In case we have a merge feature, we assume that file gives a folder location to gffs that should be merged and merge feature gives the feature to use to merge.
   if(verbose) cat ("Start to import several gffs and merge them using the feature",merge.feature,"\n")
   if(!is.element(merge.feature, features)){
     features <- c(features, merge.feature)
     cat("Added",merge.feature,"to features-Option.")
    }
   
   gffFiles <- list.files(file)
   gffFiles <- gffFiles[grepl(".gff",gffFiles)]
   gffNames <- gsub(".gff","",gffFiles)
   gffNames <- gsub(".gz","",gffNames)
   if(verbose) cat("Start to import", length(gffFiles),"files:\n", paste(gffFiles , collapse=" \n"),"\n")
   
  # Take here the feature list with names that are not used for merging
   features.small <- features[!is.element(features, merge.feature)]
   mergeThose <- c()
   
   allGFFs <- list()
   for(i in 1:length(gffFiles)){
     if(verbose) cat("Start to import: ", gffFiles[i],"\n")
     allGFFs[[i]] <- importGFF.internal(file=file.path(file,gffFiles[i]), skip=skip, nrow=nrow, use.data.table=use.data.table, level=level, features=features, num.features=num.features, print.features=print.features, verbose=verbose)
     setkeyv(allGFFs[[i]], merge.feature)
     
     takeThose <- which(is.element(colnames(allGFFs[[i]]),features.small))
     for(j in 1:length(takeThose)){
       mergeThose <- c(paste(gffNames[i],colnames(allGFFs[[i]])[takeThose[j]],sep="."))
       colnames(allGFFs[[i]])[takeThose[j]] <- mergeThose[length(mergeThose)]          
     }
     # Once the first two samples are read, merge them
     if(i==2) tmp <- merge(allGFFs[[1]], allGFFs[[2]][,unique(c(merge.feature,mergeThose)),with=FALSE], all=merge.all)     
     # Then add consecutive every turn the next one
     if(i>2){
       tmp <- merge(tmp, allGFFs[[i]][,unique(c(merge.feature,mergeThose)),with=FALSE], all=merge.all)       
     }
   }

   # Now set the columns names
   
   out <- tmp
   class(out) <- append(class(out), "gff")
 }
  out
}
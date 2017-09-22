importGTF.internal <- function(file, skip=auto, nrow=-1, use.data.table=TRUE, level="gene", features=NULL, print.features=FALSE, merge.feature=NULL, verbose=FALSE){
  gtf <- file

   if(skip=="auto"){
     if(last(strsplit(gtf,"\\.")[[1]])=="gz"){
       cat("dfg")
     }
   } else {
       con  <- file(gtf, open = "r")
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

  if(verbose) cat("Automatically detected number of rows to skip: ", skip,"\n")
  
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

    if(verbose){
      if(sum(cuffLoaded$V3==level)==0){
        availLevels <- unique(cuffLoaded$V3)
        stop("The given import level '",level,"' was not found in the gtf! Available level are: \n", paste(availLevels, collapse=" \n "))
      }
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
    if(print.features || verbose){
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

importGTF <- function(file, skip=5, nrow=-1, use.data.table=TRUE, level="gene", features=NULL, print.features=FALSE, merge.feature=NULL, class.names=NULL, verbose=FALSE){

# If no merge feature is given, we assume that only a single gtf is to be imported
   if(is.null(merge.feature)){
   out <- importGTF.internal(file=file, skip=skip, nrow=nrow, use.data.table=use.data.table, level=level, features=features, print.features=print.features, verbose=verbose)
 } else {
# In case we have a merge feature, we assume that file gives a folder location to gtfs that should be merged and merge feature gives the feature to use to merge.
   if(verbose) cat ("Start to import several gtfs and merge them using the feature",merge.feature,"\n")
   gtfFiles <- list.files(file)
   gtfFiles <- gtfFiles[grepl(".gtf",gtfFiles)]
   gtfNames <- gsub(".gtf","",gtfFiles)
   gtfNames <- gsub(".gz","",gtfNames)
   if(verbose) cat("Start to import", length(gtfFiles),"files:\n", paste(gtfFiles , collapse=" \n"),"\n")
   
  # Take here the feature list with names that are not used for merging
   features.small <- features[!is.element(features, merge.feature)]
   mergeThose <- c()
   
   allGTFs <- list()
   for(i in 1:length(gtfFiles)){
     if(verbose) cat("Start to import: ", gtfFiles[i],"\n")
     allGTFs[[i]] <- importGTF.internal(file=file.path(file,gtfFiles[i]), skip=skip, nrow=nrow, use.data.table=use.data.table, level=level, features=features, print.features=print.features, verbose=verbose)
     setkeyv(allGTFs[[i]], merge.feature)
     
     takeThose <- which(is.element(colnames(allGTFs[[i]]),features.small))
     for(j in 1:length(takeThose)){
       mergeThose <- c(paste(gtfNames[i],colnames(allGTFs[[i]])[takeThose[j]],sep="."))
       colnames(allGTFs[[i]])[takeThose[j]] <- mergeThose[length(mergeThose)]          
     }
     # Once the first two samples are read, merge them
     if(i==2) tmp <- merge(allGTFs[[1]], allGTFs[[2]][,unique(c(merge.feature,mergeThose)),with=FALSE])     
     # Then add consecutive every turn the next one
     if(i>2){
       tmp <- merge(tmp, allGTFs[[i]][,unique(c(merge.feature,mergeThose)),with=FALSE])       
     }
   }

   # Now set the columns names
   
   out <- tmp
   class(out) <- append(class(out), "gtf")
 }
  out
}
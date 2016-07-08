importPED <- function(file, verbose=TRUE){
  woEnding <- strsplit(file,".ped")[[1]][1]
  pedFile <- paste(woEnding,".ped", sep="")
  
  mapFile <- paste(woEnding,".map", sep="")
  
  if(verbose) cat("Read in the map file\n")
  map <- fread(mapFile, header=FALSE, colClasses = c("character",
                                                     "character",
                                                     "numeric",
                                                     "numeric"),
               col.names = c("Chromosome",
                             "MarkerID",
                             "GeneticDistance",
                             "Physical position")
               )

  if(verbose) cat("Read in the ped file\n")
  ped <- fread(pedFile, header=FALSE, 
               colClasses = c("numeric",
                              "character",
                              rep("numeric",4),
                              rep("character",2*nrow(map))),
               col.names = c("FamilyID",
                             "SampleID",
                             "PaternalID",
                             "MaternalID",
                             "Sex",               # (1=male; 2=female; other=unknown)
                             "Affection",          # (0=unknown; 1=unaffected; 2=affected)
                             # Genotypes (space or tab separated, 2 for each marker. 0=missing)
                             c(rbind(map$MarkerID, paste(map$MarkerID,".2",sep=""))))
                   )
  
  genotypes <- data.table(MarkerID = map$MarkerID,
                          Allele1 = rep("NA",nrow(map)),
                          Allele2 = rep("NA",nrow(map)))

  if(verbose) cat("Extract the genotype information\n")
  if(verbose) pb   <- txtProgressBar(1, nrow(map), style=3)
  for(i in 1:nrow(map)){
    all1 <- ped[[map$MarkerID[i]]]
    all2 <- ped[[paste(map$MarkerID[i],".2",sep="")]]
    alleles <- unique(c(all1,all2))
    alleles <- alleles[alleles!="0"]
    genotypes$Allele1[i] <- alleles[1]
    genotypes$Allele2[i] <- alleles[2]
    tmp <- factor(paste(all1,all2,sep=""))
    levels(tmp) <- list('0'=paste(alleles[1],alleles[1],sep=""),
                        '1'=c(paste(alleles[1],alleles[2],sep=""),
                              paste(alleles[2],alleles[1],sep="")),
                        '2'=paste(alleles[2],alleles[2],sep=""),
                        '-1'=c("00"))
   # ped[[map$MarkerID[i]]] <- as.character(tmp)
  #  ped[,paste(map$MarkerID[i],".2",sep=""):=NULL]
    if(verbose) setTxtProgressBar(pb, i)
  }
  list(ped=ped, map=map, genotypes=genotypes)
}
importVCF <- function(file, na.seq="./."){
# Necessary variable declaration for Cran checks
  V3 <- NULL
  rn <- NULL
  
 # First read in the header lines and determine the skip variable for the body
  con <- file(file) 
  open(con);
  results.list <- list();
  headerComplete <- FALSE
  headerLines <- 0
  header <- c()
  while (!headerComplete) {
    # Read in line by line
      oneLine <- readLines(con, n = 1, warn = FALSE)

    # Check if current line belongs to header
      lineStart <- substr(oneLine,1,2)
      if(lineStart=="##"){
        headerLines <- headerLines + 1
        header <- c(header, oneLine)
      } else {
        headerComplete <- TRUE
      }
    } 
  close(con)
  
  vcfBody <- fread(file, skip = headerLines, header=TRUE)
  
# Extract the map information
  map <- vcfBody[, .SD, .SDcols = c(1,3,2,4,5)]
  map[,V3:=0]
  setnames(map, c("V1", "snp.names", "V4", "allele.1", "allele.2", "V3"))
  setcolorder(map, c(1,2,6,3,4,5))
  map[[2]] <- as.character(map[[2]])
  
# If SNP names are missing, name them according to position:
  missingNames <- map[[2]]=="."
  if(sum(missingNames)>0){
    newLabels <- paste(map[[1]],map[[4]],sep=".")    
  # Test still for multi SNP per loci  
    tableNames <- table(newLabels)
    multLoci <- tableNames>1
    if(sum(multLoci)>0){
      lociOI <- tableNames[multLoci]
      for(locRun in 1:length(lociOI)){
        origLoc <- which(newLabels==names(lociOI)[locRun])
        for(indRun in 1:length(origLoc)){
          newLabels[origLoc[indRun]] <- paste(newLabels[origLoc[indRun]],indRun,sep=".")
        }
      }
    }
    map[[2]][missingNames] <- newLabels[missingNames]  
  }

  
# Extract the genotype information
  genotypes <- vcfBody[, .SD, .SDcols = -c(1:9)]
# Change them to raw format look alike, it is NOT raw!!
  for(genoRun in colnames(genotypes)){
   genotypes[get(genoRun) == na.seq, eval(genoRun) := "00"]
   genotypes[get(genoRun) == "0/0", eval(genoRun) := "01"]
   genotypes[get(genoRun) == "0/1", eval(genoRun) := "02"]
   genotypes[get(genoRun) == "1/1", eval(genoRun) := "03"]
  }
  
#  genotypes[,SNP:=map[[2]]]
  genotypes <- genotypes[, data.table(t(.SD), keep.rownames=TRUE)]  # Takes long, IMPROVE IT!!!
#  genotypes <- dcast.data.table(melt(genotypes, id.vars = "SNP"), variable ~ SNP) # This solution to above problem leads to wrong results!!!
  genotypesRN <- as.character(genotypes[[1]])
# genotypes <- genotypes[, .SD, .SDcols = -1] # TAKES LONG, IMPROVE IT!!!
  genotypes[,rn:=NULL]
  setnames(genotypes, map[[2]])
  rownames(genotypes) <- genotypesRN
  
# Then import the body
  out <- list(header=header, vcfBody, map=map, genotypes=genotypes)
  class(out) <- "vcf"
  out
}

####################################
##
## TESTING AREA FOR THE FUNCTION
##
####################################

#tmp <- importPED(file="/home/ejo138/ownCloud/Luke/Projects/Consulting/Dog-Liver/Data/NWT_151110.ped",
#                 snps="/home/ejo138/ownCloud/Luke/Projects/Consulting/Dog-Liver/Data/NWT_151110.map")
#tmp
#tmp2 <- importVCF(file="/home/ejo138/ownCloud/Luke/Projects/Consulting/Dog-Liver/Data/NWT_151110.vcf")
#tmp2

#genotDataVCF <- importVCF(file="/home/ejo138/ownCloud/R-Packages-Pages/GenomicTools/Datasets/genotypes.vcf")
#save(genotDataVCF, file="/home/ejo138/GitHub/GenomicTools/data/genotDataVCF.rda")

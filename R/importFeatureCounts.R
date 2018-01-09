importFeatureCounts <- function(file){
  tmp <- read.table(file, header=TRUE)
  expValues <- tmp[,c(1,7)]
  geneInfo <- tmp[,1:6]
  result <- list(expValues=expValues, geneInfo=geneInfo)
  class(result) <- "featureCounts"
  result
}

print.featureCounts <- function(x, ...){
  cat("$expValues \n")
  print(head(x$expValues))
  cat("...\n",nrow(x$expValues)-6,"more rows!\n")
  cat("\n")
  cat("$geneInfo")
  print(head(x$geneInfo))
  cat("...\n",nrow(x$geneInfo)-6,"morerows!\n")
}

test <- importFeatureCounts(file="/home/ejo138/ownCloud/Luke/Projects/CharacterizationOfTheRumenPapillae/counts/cow/cow_annot/UMD3.1.90/1PIntomieli2AAligned.sortedByCoord.out.bam.txt.txt")
test
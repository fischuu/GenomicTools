importFeatureCounts <- function(file){
  tmp <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
  expValues <- tmp[,c(1,7)]
  geneInfo <- tmp[,1:6]
  tmp <- read.table(paste(file,".summary",sep=""), header=TRUE, stringsAsFactors=FALSE)
  result <- list(expValues=expValues, geneInfo=geneInfo, summary=tmp)
  class(result) <- "featureCounts"
  result
}

print.featureCounts <- function(x, ...){
  cat("$expValues \n")
  print(head(x$expValues))
  cat("...\n",nrow(x$expValues)-6,"more rows!\n")
  cat("\n")
  cat("$geneInfo \n")
  print(head(x$geneInfo))
  cat("...\n",nrow(x$geneInfo)-6,"more rows!\n")
  cat("\n")
  cat("$summary \n")
  print(x$summary)
}

summary.featureCounts <- function(x, ...){
  x$summary
}

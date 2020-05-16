# This is a very early quick+dirty tool to perform that step and it needs more functions and
# performance improvements as well as input check. Further, it currently only works with two
# restriction enzymes.

restrictionDigest <- function(fa, re1="AATTC", re2="GCATG"){
  
  # Cut the reference genome at all places the fit restriction enzyme 1
  tmp <- strsplit(fa, re1)
  
  # Cut all the cutted fragments at restriction enzyme 2
  out <- list()
  for(i in 1:length(tmp)){
    out[[i]] <- strsplit(tmp[[i]],re2)
  }
  
  # Remove those fragments cutted by RE1 that did not contain a RE2 cut-site
  # Keep only those fragments that also contained a RE2 cutsite
  out.length <- list()
  cutsites <- list()
  for(i in 1:length(out)){
    out.length[[i]] <- lapply(out[[i]],length)
    cutsites[[i]] <- out[[i]][-which(unlist(out.length[[i]])==1)]
  }
  
  # Take only the first from the cusites, to ensure a RE1-----RE2 fragment, otherwise it would be RE2------RE2
  clusters <- list()
  for(i in 1:length(cutsites)){
    clusters[[i]] <- unlist(lapply(cutsites[[i]],"[",1))
  }
  
  # Now get the number of basepairs within each cluster
  clusters.bp <- list()
  for(i in 1:length(clusters)){
    clusters.bp[[i]] <- unlist(lapply(clusters[[i]],nchar))
  }
  
  out <- list(clusters=clusters,
              clusters.bp=clusters.bp)

  out  
}
# This is a very early quick+dirty tool to perform that step and it needs more functions and
# performance improvements as well as input check. Further, it currently only works with two
# restriction enzymes.

restrictionDigest <-
  function(fa,
           re1 = "AATTC",
           re2 = "GCATG",
           use.rev = TRUE) {
    if (use.rev) {
      out1 <- restrictionDigest.internal(fa = fa, re1 = re1, re2 = re2)
      out2 <-
        restrictionDigest.internal(fa = fa,
                                   re1 = revcomp(re2),
                                   re2 = revcomp(re1))
      
      out <- list(
        clusters = out1$clusters,
        clusters.bp = out1$clusters.bp,
        cutsites.re1 = out1$cutsites.re1,
        cutsites.re2 = out1$cutsites.re2,
        clusters.rev = out2$clusters,
        clusters.rev.bp = out2$clusters.bp,
        cutsites.rev.re1 = out2$cutsites.re1,
        cutsites.rev.re2 = out2$cutsites.re2
      )
    }
    else {
      out <- restrictionDigest.internal(fa = fa, re1 = re1, re2 = re2)
    }
    out
  }

restrictionDigest.internal <- function(fa, re1, re2) {
  # Cut the reference genome at all places the fit restriction enzyme 1
  tmp <- strsplit(fa, re1)
  
  cutsites.re1 <- length(tmp)
  
  # Cut all the cutted fragments at restriction enzyme 2
  out <- list()
  cutsites.re2 <- 0
  for (i in 1:length(tmp)) {
    out[[i]] <- strsplit(tmp[[i]], re2)
    cutsites.re2 <- cutsites.re2 + length(out[[i]])
  }
  
  # Remove those fragments cutted by RE1 that did not contain a RE2 cut-site
  # Keep only those fragments that also contained a RE2 cutsite
  out.length <- list()
  cutsites <- list()
  for (i in 1:length(out)) {
    out.length[[i]] <- lapply(out[[i]], length)
    cutsites[[i]] <- out[[i]][-which(unlist(out.length[[i]]) == 1)]
  }
  
  # Take only the first from the cusites, to ensure a RE1-----RE2 fragment, otherwise it would be RE2------RE2
  clusters <- list()
  for (i in 1:length(cutsites)) {
    clusters[[i]] <- unlist(lapply(cutsites[[i]], "[", 1))
  }
  
  # Now get the number of basepairs within each cluster
  clusters.bp <- list()
  for (i in 1:length(clusters)) {
    clusters.bp[[i]] <- unlist(lapply(clusters[[i]], nchar))
  }
  
  out <- list(
    clusters = clusters,
    clusters.bp = clusters.bp,
    cutsites.re1 = cutsites.re1,
    cutsites.re2 = cutsites.re2
  )
  
  out
}

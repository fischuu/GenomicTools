rearrange <- function(genoGroups,rows.gex,genoSamples){
  newOrder <- c()
  for(i in 1:length(rows.gex))
  {
    ## WARNING!!! IN CASE OF MULTIPLE MEASURES IN THE GENOTYPE, WE PICK ALWAYS THE FIRST OCCURENCE!!!!
    newOrder[i] <- which((rows.gex[i]==genoSamples)==TRUE)[1]
  }

  ifelse(is.vector(genoGroups) ,  output <- genoGroups[newOrder] , output <- genoGroups[newOrder,])
  ifelse(is.vector(genoGroups) , names(output) <- genoSamples[newOrder] , rownames(output) <- genoSamples[newOrder])
  
  output
}
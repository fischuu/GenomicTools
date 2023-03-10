# Load the data
  organism <- "ko"
  pwMap <- readLines(paste0("https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=", organism))

# Get the positions of interest
  lvl_one <- grep("<b>", pwMap)

  
output <- matrix("NA", ncol=3, nrow=length(pwMap))
current_row <- 1
for(i in 1:(length(lvl_one)-1) ){
  col_one <- pwMap[lvl_one[i]]
  col_one <- gsub("<b>", "", col_one)
  col_one <- gsub("</b>", "", col_one)
  curLevel <- 1
  for(j in (lvl_one[i]+1):(lvl_one[i+1]-1)){
    if(length(grep("<ul>", pwMap[j])) == 1) {
      curLevel <- curLevel + 1 
    } else if(length(grep("</ul>", pwMap[j])) == 1){
      curLevel <- curLevel - 1
    } else {
      if(curLevel==2){
        col_two <- trimws(pwMap[j])
      } else if(curLevel==3){
        col_three <- gsub('&nbsp;&nbsp;<a href=\"/pathway/', '', pwMap[j])
        col_three <- paste0(organism, strsplit(col_three, organism)[[1]][2])
        col_three <- gsub('\">', '_', col_three)
        col_three <- gsub('</a><br>', '', col_three)
        output[current_row,] <- c(col_one, col_two, col_three)
        current_row <- current_row + 1
      }
    }
  }
}

output <- output[1:current_row,]
output <-cbind(output, sapply(strsplit(output[,3], "_"),"[",1))
output

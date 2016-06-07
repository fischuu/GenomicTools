# This function imports the log file output from the Star Mapper

importSTARLog <- function(file){
  rawInput <- readLines(file)
  rawInput <- strsplit(rawInput," \\|\t")
  output <- cbind(sapply(rawInput,"[",1), sapply(rawInput,"[",2))
  output
}

file <- "/mnt/data2/Ruminomics/miRNA/logs/1PIntomieli2A_S1Log.final.out"
importSTARLog("/mnt/data2/Ruminomics/miRNA/logs/1PIntomieli2A_S1Log.final.out")
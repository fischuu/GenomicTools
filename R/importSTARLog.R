# This function imports the log file output from the Star Mapper

importSTARLog <- function(file){
  
  .Deprecated("GenomicTools.fileHandler::importSTARLog", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
  rawInput <- readLines(file)
  rawInput <- strsplit(rawInput," \\|\t")
  output <- cbind(sapply(rawInput,"[",1), sapply(rawInput,"[",2))
  output
}
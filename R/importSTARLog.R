# This function imports the log file output from the Star Mapper

importSTARLog <- function(file){
  rawInput <- readLines(file)
  rawInput <- strsplit(rawInput," \\|\t")
  output <- cbind(sapply(rawInput,"[",1), sapply(rawInput,"[",2))
  output
}
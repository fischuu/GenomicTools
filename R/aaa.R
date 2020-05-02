pairwiseDiffs <- function(x){
  x <- c(0,cumsum(x))
  x1 <- x[1:(length(x)-1)]
  x2 <- x[2:length(x)]
  (x1+x2)/2
}

trim.leading <- function (x)  sub("^\\s+", "", x)

chrOrder <- function(x){
  x.num <- suppressWarnings(as.numeric(x))
  nonnum.pos <- which(is.na(x.num))
  num.pos <- which(!is.na(x.num))
  
  x.num <- x[num.pos]
  x.nonnum <- x[nonnum.pos] 
  output <- c(sort(as.numeric(x.num)), sort(x.nonnum))
  
  output
}

sectoDay <- function(n){
  n <- round(n)
  
  day <- n %/% (24 * 3600) 
  n <-  n %% (24 * 3600) 
  
  hour <- n %/% 3600 
  n <- n %% 3600 
  
  minutes <- n %/% 60 
  n <- n %% 60
  
  seconds <-  n 
  
  out <- c()
  if(day > 0 ){
    out <- paste(day,"-",sprintf("%02d", hour),":",sprintf("%02d", minutes),":",sprintf("%02d", seconds),sep="")
  } else {
    out <- paste(sprintf("%02d", hour),":",sprintf("%02d", minutes),":",sprintf("%02d", seconds),sep="")
  }
  out
} 
pairwiseDiffs <- function(x){
  x <- c(0,cumsum(x))
  x1 <- x[1:(length(x)-1)]
  x2 <- x[2:length(x)]
  (x1+x2)/2
}

trim.leading <- function (x)  sub("^\\s+", "", x)
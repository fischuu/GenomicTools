# Experimental parts on time series clustering and the corresponding distance measures

calcDistances <- function(X,method="euclidean", ...){
  method <- match.arg(method, c("euclidean"))
  res <- c()
  if(method=="euclidean") res <- calcDistances.C(X)$distMat
  res
} 

# Section with new subfunction implementations

calcDistances.Splines <- function(X, spar, ...){
  NC <- ncol(X)
  NR <- nrow(X) 
# Initialize the temp list
  temp <- vector("list", NR)
  AUC <- rep(NA,NR)
  for(i in 1:NR)
  {
    temp[[i]] <- splinefun(1:100,X[1,], method="natural")
}
  #  spline2 <- splinefun(1:100,X[2,], method="natural")
    
  #  plot(1:100,X[1,])
#    curve(spline2(x), from=1, to=10, col=1, n=1001)
#    curve(spline1(x), add=TRUE, col=2, n=1001)
#    curve(spline2(x), add=TRUE, col=3, n=1001)
#    abline(0,0)
#    points(X[1,])
#    points(X[2,])
#    
#    integrate(temp[[1]], lower=1, upper=100, subdivisions=1000)$value
#  }
#  temp
}

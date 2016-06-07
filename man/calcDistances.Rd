\name{calcDistances}
\alias{calcDistances}

\title{
  Calculate Distances Between Time Series
}

\description{
  This is an experimental function to calculate distances between time series.
}
\usage{
  calcDistances(X, method="euclidean", ...)
}

\arguments{
  \item{X}{Matrix with time series values.}
  \item{method}{Method  of choice to calculate the distance.}
  \item{...}{Additional values.}
}

\details{
 This is an early attempt to implement a distance measure for time series clustering. It has not been used in
 productive code yet and hence should be used carefully.

 Currently there is only the basic euclidean method available the default is \code{method="euclidean"} 
 
}
\value{
 A matrix with pairwise distances.
}

\author{
  Daniel Fischer
}

\keyword{
  methods
}
 

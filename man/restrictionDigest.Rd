\name{restrictionDigest}
\alias{restrictionDigest}

\title{
  Simulate DNA fragmentation with restriction enzymes.
}

\description{
  This function determines the expected fragments from a DNA fragmentation with restriction enzymes.
}
\usage{
  restrictionDigest(fa, re1="AATTC", re2="GCATG", use.rev=TRUE)
}

\arguments{
  \item{fa}{A reference genome in fasta format}
  \item{re1}{The first restriction enzyme sequence}
  \item{re2}{The second restriction enzyme sequence}
  \item{use.rev}{Logical, should the reverse complement be used}
}

\details{
 This is a very early quick+dirty tool to perform that step and it needs more functions and
 performance improvements as well as input check. Further, it currently only works with two
restriction enzymes.
}

\value{
 A list containing the clusters (\code{clusters}) and their lengths (code{cluster.bp})
}


\author{
  Daniel Fischer
}

\keyword{
  multivariate
}
 

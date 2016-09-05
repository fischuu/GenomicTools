\name{importVCF}
\alias{importVCF}

\title{
  Import a vcf file.
}

\description{
  This function imports a vcf filepair.
}
\usage{
  importVCF(file, na.seq="./.")
}

\arguments{
  \item{file}{The vcf filename.}
  \item{na.seq}{Sequence, how missing genotypes are encoded.}
}

\details{
 This function imports a vcf file for further analysis. 
}
\value{
 A \code{vcf} object
}


\author{
  Daniel Fischer
}

\keyword{
  multivariate
}
 

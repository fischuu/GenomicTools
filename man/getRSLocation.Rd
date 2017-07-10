\name{getRSLocation}
\alias{getRSLocation}

\title{
  Get the location of an SNP.
}

\description{
  This function gets the exact location of a SNPs, based on a specific database.
}
\usage{
  getRSLocation(rs, species="Gallus_gallus", url="jul2016.archive.ensembl.org")
}

\arguments{
  \item{rs}{SNP rs number}
  \item{species}{Species identifier}
  \item{url}{Path of the ensembl database to use}
}

\details{
  This function can download from any specifi ensembl database the location of an SNP (as given as RS number) for a specific species.
}
\value{
  A list with Chr and Location (BP) of an SNP.
}

\examples{
  getRSLocation("rs16033432")
}

\author{
  Daniel Fischer
}

 

\name{importPED}
\alias{importPED}

\title{
  Import a Ped/Map filepair.
}

\description{
  This function imports a ped/map filepair.
}
\usage{
  importPED(file, n, snps, which, split = "\t| +", sep = ".",
            na.strings = "0", lex.order = FALSE, verbose=TRUE)
}

\arguments{
  \item{file}{The ped filename.}
  \item{n}{Number of observations to import.}
  \item{snps}{The map file name.}
  \item{which}{A set of SNPs to import.}
  \item{split}{Column separator in ped/map filepair.}
  \item{sep}{Allele separator.}
  \item{na.strings}{Label for missing data.}
  \item{lex.order}{Logical, shall the data be ordered}
  \item{verbose}{Logical, gives verbose feedback}
}

\details{
 This function is to a large extend taken from \code{snpStat::read.pedmap}, but here is internally the \code{data.table::fread} function used
 that resulted in much faster file processing. 
}
\value{
 A \code{PedMap} object
}


\author{
  Daniel Fischer
}

\keyword{
  multivariate
}
 

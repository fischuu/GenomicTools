\name{gtfToBed}
\alias{gtfToBed}

\title{
  Extract the Chromosomal Information Required in bed Format from an imported gtf table.
}

\description{
  This function creates a matrix of gene annotations in bed format, based on the information given in an importet gtf table.
}
\usage{
  gtfToBed(gtf)
}

\arguments{
  \item{gtf}{An imported gtf table.}
}

\details{
 After a gtf file was imported with \code{importGTF} this function transforms it into a data.frame with bed formatting.
}
\value{
 A data.frame in bed format having the four columns \code{Chr}, \code{Start}, \code{Stop} and \code{Gene}
}


\author{
  Daniel Fischer
}

\examples{
\dontrun{
annotTrack <- read.table(file="Homo_sapiens.GRCh37.70.gtf",sep="\t")
annotBed <- gtfToBed(annotTrack)
}
}

\keyword{
  multivariate
}
 

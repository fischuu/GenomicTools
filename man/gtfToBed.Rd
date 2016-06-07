\name{gtfToBed}
\alias{gtfToBed}

\title{
  Extract the Chromosomal Information Required in bed Format from an imported gtf table.
}

\description{
  This function creates a matrix of gene annotations in bed format, based on the information given in an importet gtf table.
}
\usage{
  gtfToBed(gtf, output="min")
}

\arguments{
  \item{gtf}{An imported gtf table.}
  \item{output}{Option to format the output, see details.}
}

\details{
  Currently the function supports only gtf files for human organisms. If applied to other organisms the code has to be changed in such a way that the Chromosome names are adjusted.
  The output can be formated with the output option, where \code{output="min"} results in one row per gene and \code{output="full"} keeps the original format in the gtf file.
}
\value{
 A data.frame in bed format having the four columns \code{Chr}, \code{Start}, \code{Stop} and \code{Name}
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
 

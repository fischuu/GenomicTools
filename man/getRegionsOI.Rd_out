\name{getRegionsOI}
\alias{getRegionsOI}

\title{
  Filter Annotation File
}

\description{
  Filter out unimportant regions of an annotation file in bed format. 
}
\usage{
   getRegionsOI(annot,regOI)
}

\arguments{
  \item{annot}{An annotation file in bed format.}
  \item{regOI}{A dataframe specifying the interesting location.}
}

\details{
 When considering trans-eQTLs all combinations between gene annotations and genotype data are tested. This often forces the user
 to filter the annotations to a certain region. Having an annotation file in bed format with the columns \code{Name},\code{Chr}, \code{Start}
 and \code{Stop} and a data frame specifying the different locations that should remain, this function removes all annotations outside these 
 given areas. THe column names of the data frame given to \code{regOI} are \code{Chr}, \code{Start} and \code{End}.
}
\value{
 A subset of \code{annot}.
}

\author{
  Daniel Fischer
}

\keyword{
  methods
}
 

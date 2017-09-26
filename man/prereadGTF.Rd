\name{prereadGTF}
\alias{prereadGTF}

\title{
  Pre-read a GTF file.
}

\description{
  This function imports a gtf file and prints features of it for importing it.
}
\usage{
  prereadGTF(file, nrow=1000, skip="auto")
}

\arguments{
  \item{file}{The gtf filename.}
  \item{nrow}{Number of rows to be read}
  \item{skip}{Number of rows to skip}
}

\details{
 This function reads in a gtf file and prints its features for the import step. 
 
 By default this function only imports the first 1000 rows, in case all rows should be imported set \code{nrow=-1}.
 
 The number to skip in the beginning can be adjusted by the \code{skip} option. The default is here \code{auto} so that
 the function can identify the correct amount of header rows. Hence, this option should be changed only, if there is a
 good reason.
 
}
\value{
 A summary of the gtf features
}


\author{
  Daniel Fischer
}

\seealso{
\code{importGTF}
}

\keyword{
  multivariate
}
 

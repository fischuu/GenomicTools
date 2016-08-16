\name{mdr}
\alias{mdr}

\title{
  Perform a MDR.
}

\description{
  This function performs a Multifactor Dimensionality Reduction (MDR).
}
\usage{
  mdr(X, status, fold=2, t=NULL, top=3, NAasValues=TRUE, 
      fix=NULL, verbose=FALSE)
}

\arguments{
  \item{X}{Matrix with genotype information, see details.}
  \item{status}{Vector with group information of individuals in \code{X}.}
  \item{fold}{Maximum dimension of used contingency tables, see details.}
  \item{t}{Threshold for high/low risk.}
  \item{top}{Length of each top list.}
  \item{NAasValues}{How shall NAs be treated.}
  \item{fix}{Column number of the SNP to be fixed.}
  \item{verbose}{Logical, if detailed status messages shall be given.}
}

\details{
 The matrix \code{X} contains the genotype information, each column corresponds to a SNP, each row to an individuum. SNPs need to be coded as 0,1,2.
 In case the matrix \code{X} is not given in 0,1,2 format the function \code{recodeData} recodes the data into the required 0,1,2 format.

 The \code{status} vector is as long as \code{X} has individuals/rows and specifies the group labels for each individual. Healthy individual need to be
 encoded as 0 and cases as 1. If the labeling is different the smaller values are used as controls and the larger one as cases. 

 The \code{fold} option specifies up to which dimension the contingency tables should be used. The current maximum is four.

 The \code{t} option gives the threshold for the classification between high and low risk classes. The default
 is the ratio of the groups sizes.
 
 With the \code{top} option the amount of results are set. 
 
 Setting the \code{fix} option to a column number, forces the mdr function to include that particular SNP into the result.
}
\value{
 An object of class \code{mdr}.
}


\references{
  \cite{Moore JH, Gilbert JC, Tsai CT, Chiang FT, Holden T, Barney N, White BC. (2006): A flexible computational
   framework for detecting, characterizing, and interpreting statistical patterns of epistasis in genetic studies
   of human disease susceptibility. J Theor Biol.2006 Jul 21;241(2):252-61.}

}

\author{
  Daniel Fischer
}

\examples{
data(mdrExample)
mdrSNP <- mdrExample[,1:20]
fit.mdr <- mdr(mdrSNP, mdrExample$Class, fold=4, top=5)
fit.mdr
fit.mdr <- mdr(mdrSNP, mdrExample$Class)
fit.mdr
}

\keyword{
  multivariate
}
 

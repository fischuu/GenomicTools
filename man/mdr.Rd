\name{mdr}
\alias{mdr}

\title{
  Perform a MDR.
}

\description{
  This function performs a Multifactor Dimension Reduction (MDR).
}
\usage{
  mdr(X,status,fold=2,t=NULL,cv=0,cvp=0.75,top=20,NAasValues=TRUE,fix=NULL)
}

\arguments{
  \item{X}{Matrix with genotype information, see details.}
  \item{status}{Vector with group information of individuals in \code{X}.}
  \item{fold}{Maximum dimension of used contingency tables, see details.}
  \item{t}{Threshold for high/low risk.}
  \item{cv}{Amount of cross validation runs.}
  \item{cvp}{Ratio of cross-validation sample.}
  \item{top}{Length of each top list.}
  \item{NAasValues}{How shall NAs be treated.}
  \item{fix}{Shall one genotype be fixed.}

}

\details{
 The matrix \code{X} contains the genotype information or the filename of a ped/map filepair. If a ped/map
 filename is given the status information from this pair is taken and no further \code{status} object has
 to be given. In case \code{status} is given as well, we will take this information.

 In case the matrix \code{X} is not given in 0,1,2 format the function \code{recodeData} recodes the data 
 into the required 0,1,2 format.

 The \code{status} vector is as long as \code{X} has individuals and specifies the group labels for each
 individual. Healthy individual shall be encoded as 0 and cases as 1. If the labeling is different we take 
 the smaller values as controls and the larger one as cases. 

 The \code{fold} option specifies up to which dimension the contingency tables should be used. The current 
 maximum is four.

 The \code{t} option gives the threshold for the classification between high and low risk classes. The default
 is the ratio of the groups sizes.
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
# The datasets are not yet available. As soon as they will get published they will
# be also added to the package

\dontrun{
# Read in the genotype data
genotData <- read.table("MDR_data.txt",header=T)

# Extract the status information
status <- genotData[["Class"]]
genotData <- as.matrix(genotData[,-which((colnames(genotData)=="Class")==TRUE)])

# Bring the data into 0,1,2 format
temp <- recodeData(genotData)

# Perform the MDR
res <- mdr(X=temp,status=status,fold=3,top=20)
}
}

\keyword{
  multivariate
}
 

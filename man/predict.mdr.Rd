\name{predict.mdr}
\docType{methods}
\alias{predict.mdr}

\title{Calculates MDR Ensemble Classifier}
\description{
  This is an mdr ensemble classifier. 
}
\usage{
 \method{predict}{mdr}(object, data=NULL, status=NULL, fold=NULL, ...)
}

\arguments{
  \item{object}{Object of class \code{mdr}.}
  \item{data}{The new data object.}
  \item{status}{Optional, used for 2x2 classification table.}
  \item{fold}{Considered dimension of the model.}
  \item{...}{Additional parameters}
}

\details{
Given an \code{mdr} object this function takes the top list for the (highest) fold group and uses it as
an ensemble classifier for the new data given in the \code{data} argument.
}

\value{
A vector, giving for each subject from the new data object a classification.
}

\examples{
# indices <- 1:nrow(genotData)
# trainSet <- sample(indices,100)
# testSet <- indices[-trainSet]
# temp <- recodeData(genotData)
# res <- mdr(X=temp[trainSet,], status=status[trainSet], fold=3, top=20)

# trainRes <- predict(res, data=temp[trainSet,])
# testRes <- predict(res, data=temp[testSet,])
}

\author{Daniel Fischer}

\keyword{methods}

\name{mdrEnsemble}
\alias{mdrEnsemble}

\title{
  Ensemble Classification Using an Existing MDR Object
}

\description{
  This function takes an exisitng MDR object and a set of new observations and classifies them according to the best n serts of MDR classifiers.
}
\usage{
  mdrEnsemble(mdr, data=NULL, new.status=NULL, fold=NULL)
}

\arguments{
  \item{mdr}{The MDR object.}
  \item{data}{The new data.}
  \item{new.status}{Optional, the true status of the new data.}
  \item{fold}{Level of SNP-SNP interaction.}
}

\details{
 This function takes an exisitng MDR object and a set of new observations and classifies them then using the top results of it, as it was defined
 in the previous MDR run.  
}
\value{
 A classification of the new observations
}


\author{
  Daniel Fischer
}

\keyword{
  multivariate
}
 

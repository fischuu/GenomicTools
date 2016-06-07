\name{visTrans}
\alias{visTrans}

\title{
  Visualize trans-eQTL Results
}

\description{
  Plot function for the visualization of trans-eQTL results.
}
\usage{
  visTrans(snpGene, geneAnnot)
}

\arguments{
  \item{snpGene}{A dataframe indicating the snp-gene association, see details.}
  \item{geneAnnot}{An annotation track.}
}

\details{
 Typically are trans-eQTL difficult to visualize. One possible option is this plot. It takes 
 as an input a dataframe in bed format. This one is typically provided from an \code{eQTL} run
 with set parameter \code{sig} in the \code{bed} list. The dataframe indicates then all significant 
 associations between genes and SNPs and these associations are then connected within the Figure
 with arches. 
 
 In addition it is possible to plot an annotation track, by specifying it in \code{geneAnnot}. This 
 track is required to be in bed format and the column names are  \code{Name}, \code{Chr}, \code{Start}
 and \code{Stop}. In case a standard Ensemble gtf file is used for that the function \code{gtfToBed}
 provides the correct input for this option.
}
\value{
 A figure.
}

\author{
  Daniel Fischer
}

\examples{
  # This is just simulated data and hence we cannot see results as with real data. 
  # An example for real data is shown in the vignette and as soon as this data is
  # freely available this example will be updated.
  
  # See also the vignette for a more detailed example.
 
  data(Xgene) 
  data(genotData)
  data(annotTrack)
\dontrun{
  annotBed <- gtfToBed(annotTrack)
  
  lm.myEQTLtrans <- eQTL(gex=Xgene, geno=genotData, xAnnot=annotBed, 
                         method="LM", windowSize=NULL, sig=0.01)

  snpGeneInfo <- lm.myEQTLtrans$bed
# Plot the visualization
  visTrans(snpGene=snpGeneInfo,annotBed)
}
}

\keyword{
  hplot
}
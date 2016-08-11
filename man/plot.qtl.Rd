\name{plot.qtl}
\docType{methods}
\alias{plot.qtl}
\alias{plot-method}
\alias{plot,qtl-method}

\title{Plot an qtl Object}
\description{
 The function offers informative plots for an \code{qtl} object.
}
\usage{
 \method{plot}{qtl}(x, which=NULL, sig=0.01, verbose=TRUE,
                    log=FALSE, genome=NULL, pch=1, ...)
 }

\arguments{
  \item{x}{Object of class \code{qtl}.}
  \item{which}{Specifies for which phenotypes the plot should be created.}
  \item{sig}{Chosen significance level.}
  \item{verbose}{Logical, extended feedback of the function.}
  \item{log}{Logical, y-axis scale is log(base=10)-scaled.}
  \item{genome}{Data frame with Chromosome names and lengths.}
  \item{...}{Additional plotting parameters.}
}

\details{
 This function plots the test results of an \code{eqtl} object. Typically is the tested gene in the center and the p-values of associated SNPs are visualized. 
 Monomorphic SNPs and those that were missing are separately plotted. Test results that are smaller than the value given to \code{sig} are marked in red. 
 The y.axis can be switched to log10 scale by setting the logical parameter \code{log=TRUE} in that case are bars instead of dots plotted. If the y-axis is on log-scale
 it is also possible to give a second \code{eqtl} object to the function and plot the test results for both.
 
 The annotation feature is currently under development and only available in limited form. 
}

\author{Daniel Fischer}

\examples{
# Perform eQTL (single location, one gene):
  data(Xgene) 
  data(genotData)
  data(annotTrack)
\dontrun{
  annotBed110 <- gtfToBed(annotTrack[1:10,])

  lm.myEQTL <- eQTL(gex=Xgene,geno=genotData, xAnnot=annotBed110,method="LM")
  
  plot(lm.myEQTL)
  
  dir.myEQTL <- eQTL(gex=Xgene,geno=genotData, xAnnot=annotBed110,method="directional")
  plot(lm.myEQTL, x2=dir.myEQTL, log=TRUE, sig=2)
}
}

\keyword{methods}
\keyword{hplot}

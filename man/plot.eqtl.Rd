\name{plot.eqtl}
\docType{methods}
\alias{plot.eqtl}
\alias{plot-method}
\alias{plot,eqtl-method}

\title{Plot an eqtl Object}
\description{
 The function offers informative plots for an \code{eqtl} object.
}
\usage{
 \method{plot}{eqtl}(x, file = NULL, which = NULL, sig = 0.01, verbose = TRUE, centered = TRUE,
                     log = FALSE, x2 = NULL, annot = NULL, track = NULL,
                     trackAnnot = FALSE, trackOrder = NULL, mc.cores = 1, ...)
}

\arguments{
  \item{x}{Object of class \code{eqtl}.}
  \item{file}{Store set of graphics under that file name.}
  \item{which}{Specifies for which genes should the plot be created.}
  \item{sig}{Chosen significance level.}
  \item{verbose}{Logical, extended feedback of the function.}
  \item{centered}{Logical, plot should be centered around center gene.}
  \item{log}{Logical, y-axis scale is log(base=10)-scaled.}
  \item{x2}{Comparison values of a second \code{eqtl} object, see details.}
  \item{annot}{Logical, plot annotation track.}
  \item{track}{Gene annotations in bed format.}
  \item{trackAnnot}{Gene annotations in bed format.}
  \item{trackOrder}{Logical, shall the annotation track be ordered.}
  \item{mc.cores}{Amount of cores for parallel computing.}
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

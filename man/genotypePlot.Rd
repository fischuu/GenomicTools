\name{genotypePlot}
\alias{genotypePlot}

\title{
  Expression Boxplots
}

\description{
  Expression values are grouped according to genotype groups and then visualized with boxplots.
}
\usage{
  genotypePlot(snp,gene,eqtl,geneAnnot=NULL, ylab=NULL, xlab=NULL, mainlab=FALSE)
}

\arguments{
  \item{snp}{\code{String}, specifies the genotype name.}
  \item{gene}{\code{String}, specifies the gene name.}
  \item{eqtl}{An \code{eqtl} object.}
  \item{gex}{Matrix or vector with expression values}
  \item{geno}{Genotype data}
  \item{ylab}{Optional x-axis label.}
  \item{xlab}{optional y-axis label}
  \item{mainlab}{\code{Logical}, shall main title be plotted}
}

\details{
 This function plots the expression values of the genotype groups of a certain SNP that can be given in the
 \code{snp} option. The expression values are specified in the \code{gene} option. 
 
 The \code{eqtl} object is the output of an \code{eQTL} run and carries the required genotype information.  
}
\value{
 A Figure.
}

\author{
  Daniel Fischer
}

\examples{

# Make the example data available
  data(Xgene) 
  data(genotData)
  data(annotTrack)

# We need to have the gene annotation in bed format (Please notice the change to the
# official convention, this is on high priority of the ToDo list of the package to change
# this.)
\dontrun{
  annotBed <- gtfToBed(annotTrack)
  
# Perform a basic cis-eQTL with the minimum required input linear model:
  lm.myEQTL <- eQTL(gex=Xgene,geno=genotData, xAnnot=annotBed,method="LM",windowSize=1)

# Plot the genotypes
  genotypePlot(snp="SNP377", gene="MYBPC1", eqtl=lm.myEQTL, ylab="Expression values", 
               xlab="Genotypes")
}
}

\keyword{
  hplot
}
 

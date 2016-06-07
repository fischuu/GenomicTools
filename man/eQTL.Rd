\name{eQTL}
\alias{eQTL}

\title{
  Perform an eQTL Analysis
}

\description{
  This function performs an eQTL analysis.
}
\usage{
  eQTL(gex, geno, xAnnot = NULL, xSamples = NULL, genoSamples = NULL,
       windowSize = 0.5, method = "LM", mc = 1, sig = NULL, which = NULL,
       nper = 2000, verbose = TRUE)
}

\arguments{
  \item{gex}{Matrix or Vector with expression values.}
  \item{geno}{Genotype data.}
  \item{xAnnot}{Location annotations for the expression values.}
  \item{xSamples}{Sample names for the expression values, see details (optional).}
  \item{genoSamples}{Sample names for the genotype values, see details (optional).}
  \item{windowSize}{Size of the window around the center gene, see details.}
  \item{method}{Method of choice for the eQTL, see details.}
  \item{mc}{Amount of cores for parallel computing.}
  \item{sig}{Significance level for the eQTL testing, see details.}
  \item{which}{Names of genes for that the eQTL should be performed.}
  \item{nper}{Sets the amount of permutations, if permuation tests are used.}
  \item{verbose}{Logical, if the method should report intermediate results.}
}

\details{
This function performs an eQTL analysis and offers different types of tests. The type of test
can be specified with the \code{method} option and possible options are \code{"LM"} and \code{"directional"}. 
The option \code{"LM"} fits for each SNP within a predefined window of size \code{windowSize} (in MB) around a gene
a linear model for the genotype information and the corresponding gene expression. The null hypothesis
for each test is then that the slope is equal to zero and the alternative is that it is not zero.

The \code{"directional"} option applies a new directional test based on probabilistic indices for triples as described 
in Fischer, Oja, et al. (2014). Being \eqn{\mathbf{x}_0=(x_{01},x_{02},\dots,x_{0N_0})', \mathbf{x}_1=(x_{11},x_{12},\dots,x_{1N_1})'}
and \eqn{\mathbf{x}_2=(x_{21},x_{22},\dots,x_{2N_2})'} the expression values that are linked to the three genotype
groups \eqn{0,1} and \eqn{2} with underlying distributions \eqn{F_0, F_1} and \eqn{F_2}.  We first calculate the probabilistic
indices \eqn{P_{0,1,2} = \frac{1}{N_0 N_1 N_2} \sum_i \sum_j \sum_k I(x_{0i} < x_{1j} < x_{2k})}
and \eqn{P_{2,1,0} = \frac{1}{N_0 N_1 N_2} \sum_i \sum_j \sum_k I(x_{2i} < x_{1j} < x_{0k})}. These are the probabilities that the expression
values of the three groups follow a certain order what we would expect for possible eQTLs. The null hypothesis that we have then
in mind is that the expression values from these three groups have the same distribution \eqn{H_0: F_0 = F_1 = F_2} and the
two alternatives are that the distributions have a certain stochastical order \eqn{H_1: F_0 < F_1 < F_2} and \eqn{H_2: F_2 < F_1 < F_0}.

 The test is applied for the two probabilistic indices \eqn{P_{0,1,2}} and
\eqn{P_{2,1,0}} and combines the two resulting p-values \eqn{p_{012}=p_1} and \eqn{p_{210}=p_2} from previous tests then as overall p-value
\eqn{\min(2 \min(p_1 , p_2 ), 1)}. In the two-group case (this means only two different genotypes are present for a certain SNP) a two-sided
Wilcoxon rank-sum test is applied.

The gene expressions used in the eQTL are specified in \code{gex}. If several genes should be tested, then \code{gex} is a matrix and each
column refers to a gene and each row to an individuum. The column names of this matrix must match then with the 
names used in the annotation object \code{xAnnot}. Sample names can either be given as row names in the matrix or as separate
vector in \code{xSamples}. If only gene expressions of one gene should be tested then \code{gex} can be a vector.

The genotype information is provided in the \code{geno} object. Here one can i.a. specify the file
name of a ped/map file pair. The function then imports the genotype information using the
\code{SnpStats} package. In case the genotype information has been imported already earlier using
\code{SnpStats::read.pedfile()} the resulting \code{SnpMatrix} can also be given as a parameter for \code{geno}.

The \code{xAnnot} object carries the annotation information for the gene expressions. In case of multiple locations per gene 
it is of type list and each list item stores the information for one gene in form of a \code{data.frame} in bed format. This \code{data.frame}
has then the three columns \code{Chr}, \code{Start}, \code{End} and each row refers to one matching chromosomal postion of the underlying gene.
Especially when probes of ssRNAs are considered the chromosomal positions of a probe are not necessarily unique. The names
of the list \code{xAnnot} are the names of the genes and they have to match with the column names of \code{gex}. However, the order
does not have to be the same, and \code{xAnnot} can include more annotations of genes than given in \code{gex}. The function finds
and uses then the union between the column names of \code{gex} and the list entries of \code{xAnnot}.
Alternative \code{xAnnot} can also be a data frame if unique locations are considered. In that case \code{xAnnot} has
to be a data frame with the four columns \code{Gene}, \code{Chr}, \code{Start}, \code{End}.

The option \code{genoSamples} is used in case that the sample names in the ped/map file (or SnpMatrix)
do not match with \code{rownames(gex)} given in the expression matrix. The vector \code{genoSamples} is as
long as the \code{geno} object has samples, but gives then for each row in \code{geno} the corresponding
name in the \code{gex} object. The function finds then also the smallest union between the two data objects.
If there are repeated measurements per individual for the genotypes we take by default only the first
appearance in the data and neglect all successive values. Currently this cannot be changed. In case
this behavior is not desired, the user has to remove the corresponding rows from \code{geno} before starting the 
calculation.

If the code is executed on a Linux OS the user can specify with the \code{mc} option the amount of CPU cores used for the
calculation. 

If the \code{sig} option is set to a certain significance level, then the method only reports those SNPs that are tested to be
significant. This can reduce the required memory drastically, especially in the case of trans-eQTL. 

The method tests for trans-eQTLs (all combinations of SNPs and genes) if the \code{windowSize} is set to \code{0} or \code{NULL}.
Be aware that this might lead to long lasting calculations. 

Note: The directional test currently supports only exact p-values based on permutation tests, but asymptotic implementations are developed and 
will be soon available also.

}
\value{
 A list of class \code{eqtl} containing the values 
 \item{gex}{The \code{gex} object from the function call.}
 \item{geno}{The \code{geno} object from the function call.}
 \item{xAnnot}{The \code{xAnnot} object from the function call.}
 \item{genoSamples}{The \code{genoSamples} object from the function call.}
 \item{windowSize}{The \code{windowSize} object from the function call.}

 and an incapsulated list \code{eqtl} where each list item is a tested gene location and contains the items

 \item{ProbeLoc}{Used position of that gene. (Only different from 1 if multiple locations are considered.)}
 \item{TestedSNP}{Details about the considered SNPs.}
 \item{p.values}{P values of the test.}
 \item{GeneInfo}{Details about the center gene.}
}

\references{
  \cite{Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2014): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistics, 3, pages 672-692, 
           doi: 10.1111/sjos.12055.}

  \cite{Fischer, D., Oja, H. (2013): Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT, submitted article.}
}


\author{
  Daniel Fischer
}

\examples{

# Please, see also the package vignette for a more descriptive example section on this.

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
}
}

\keyword{
  methods
}
 

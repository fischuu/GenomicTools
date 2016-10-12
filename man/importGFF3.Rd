\name{importGFF3}
\alias{importGFF3}

\title{
  Import a GFF3 File
}

\description{
  This function imports a gff3 file, as provided by NCBI.
}
\usage{
  importGFF3(gff, chromosomes)
}

\arguments{
  \item{gff}{File name of the gff3 file.}
  \item{chromosomes}{Matrix with chromosome/accession number pair.}
}

\details{
  This function imports a standard gff3 file as provided by NCBI. The gff3 files on NCBI do not have a chromosome identifier,
  but they use the accession ID of the chromosome instead. For that reason, on can provide a matrix (as 'chr_accessions...' on
  the NCBI ftp server). 
}
\value{
  A data table.
}

\author{
  Daniel Fischer
}

 

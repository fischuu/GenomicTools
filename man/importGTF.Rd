\name{importGTF}
\alias{importGTF}

\title{
  Import a GTF File
}

\description{
  This function imports a gtf file.
}
\usage{
  importGTF(file, skip = auto, nrow = -1, use.data.table = TRUE, level = "gene",
            features = NULL, num.features=c("FPKM", "TPM"), print.features=FALSE,
            merge.feature=NULL, verbose=TRUE)
}

\arguments{
  \item{file}{File name of the gtf file.}
  \item{skip}{Rows to skip from the top (usualy the number of header rows.}
  \item{nrow}{Total amount of rows read.}
  \item{use.data.table}{Logical, shall the data.table package be used for importing the GTF file.}
  \item{level}{Filter for the GTF Entry.}
  \item{features}{Feature names to extract from column #9.}
  \item{num.features}{Feature names that contain numeric values.}
  \item{print.features}{logical, shall a list of availbe features be printed.}
  \item{merge.feature}{The name of the feature to merge the data.}
  \item{class.names}{Vecotr with labels for the merged features.}
  \item{verbose}{Logical, for 'talk-active' output}
}

\details{
  This function imports a gtf file. The features names to be imported are defined in \code{features}, several features are then
  provided as vector. A list of available feature can beprinted, by setting \code{print.features=TRUE}.
  
  The \code{skip} option allows to skip a given number of rows, the default is, however, \code{auto}. In that case, all rows that
  start with the \code{#} symbol are skipped.
  
  In case a set of expression values given in gtf format should be imported and to be merged into a single data table, the feature
  that should be used for merging can be provided to the \code{merge.feature} option. In that case the function expects a folder
  in \code{file} and it will import all gtfs located in that folder and merges them according to the \code{merge.feature} option.
  With the option \code{class.names} a vector of prefixes for the merged features can be provided. If this is kept empty, then the 
  filenames of the gtf will be used instead (without gtf extension).
  
  By default the function imprts all features in column 9 as string character. However, for common labels (FPKM and TPM) the class
  type is set automatically to numeric. Additional numerical feature names can be defined with the \code{num.feature} option.
  
}
\value{
  A data table of class \code{gtf} containing the columns of the gtf file, including the splitted last column.
}

\author{
  Daniel Fischer
}

 

\name{importGTF}
\alias{importGTF}

\title{
  Import a GTF File
}

\description{
  This function imports a gtf file.
}
\usage{
  importGTF(gtf, skip = 0, nrow = -1, use.data.table = FALSE)
}

\arguments{
  \item{gtf}{File name of the gtf file.}
  \item{skip}{Rows to skip from the top.}
  \item{nrow}{Total amount of rows read.}
  \item{use.data.table}{Logical, shall the data.table package be used for importing the GTF file.}
}

\details{
  This function imports a gtf file and splits the column 9 which is usually tricky to handle as the order of the variables is not always
  the same.
  
  This function is currently improvement and its use might be unefficient.
}
\value{
  A data frame containing the columns of the gtf file, including the splitted last column.
}

\author{
  Daniel Fischer
}

 

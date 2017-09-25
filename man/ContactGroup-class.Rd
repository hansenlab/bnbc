\name{ContactGroup-class}
\Rdversion{1.1}
\docType{class}
\alias{ContactGroup-class}
\alias{[,ContactGroup,ANY,ANY,ANY-method}
\alias{colData,ContactGroup-method}
\alias{colData<-,ContactGroup,DataFrame-method}
\alias{dim,ContactGroup-method}
\alias{rowData,ContactGroup-method}
\alias{rowData<-,ContactGroup-method}
\alias{show,ContactGroup-method}
\alias{ContactGroup}
\alias{contacts}
\alias{contacts<-}
\alias{librarySize}
\alias{logCPM}

\title{Class \code{"ContactGroup"}}
\description{
  The ContactGroup class represents a collection of contact matrices
  which are observations from different samples on the same set of
  genomic loci.
}
\usage{
ContactGroup(rowData, contacts, colData)
}
\arguments{
  \item{rowData}{Object of class \code{GenomicRanges} equal
    in length to the number of rows/columns in contact matrices. }
  \item{contacts}{Object of class \code{list} that contains
    all contact matrices.}
  \item{colData}{Object of class \code{DataFrame} containing
    sample-level information.}
}

\section{Methods}{
  In the code snippets below, \code{x} is a ContactGroup object.
  \describe{
    \item{[}{\code{signature(x = "ContactGroup", i = "ANY", j = "ANY",
	drop = "ANY")}: Allows for subsetting the contact matrices
      through use of \code{i} or of samples through \code{j}.}
    \item{colData}{\code{signature(x = "ContactGroup")}: Get
      sample-level information about samples in \code{x}}
    \item{colData<-}{\code{signature(x = "ContactGroup", value =
	"DataFrame")}: Set sample-level information about samples in \code{x}.
      \code{value} is expected to be a \code{DataFrame} object.} 
    \item{dim}{\code{signature(x = "ContactGroup")}: Obtain the
      dimensions of a ContactGroup.  Returns 2 values: one 
      representing the number of bins in the contact matrices and
      another representing the number of samples.} 
    \item{rowData}{\code{signature(x = "ContactGroup")}: Get the
      GenomicRanges object describing the loci in the ContactGroup.
      \code{value} is expected to be a GenomicRanges object.}
    \item{rowData<-}{\code{signature(x = "ContactGroup")}: Set the
      GenomicRanges object describing the bins in the ContactGroup.
      \code{value} is expected to be a GenomicRanges object.}
    \item{show}{\code{signature(object = "ContactGroup")}: Method to
      display summary information about a ContactGroup: the number of
      bins, the width of the bins and the number of samples.} 
  }
}
\section{Utilities}{
  \describe{
    \item{}{\code{contacts(x)}, \code{contacts(x) <- value}: Method to
      extract the list of contact matrices from a
      ContactGroup. \code{value} is expected to be a list object.}
  }
}
\examples{
data(chr22)

chr22[1,]
chr22[,1]

cd <- colData(chr22)
colData(chr22) <- cd

gr <- rowData(chr22)
rowData(chr22) <- gr

chr22

cl <- contacts(chr22)
contacts(chr22) <- cl
}
\keyword{classes}
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
\alias{distanceIdx}
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
      drop = "ANY")}: 
      Allows for subsetting the contact matrices
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
    \item{librarySize}{\code{signature(x = "ContactGroup")}: Method to
      compute the library size of each contact matrix in \code{x}.
      Library size is defined to be the sum of the upper triangle of a
      contact matrix.}
    \item{logCPM}{\code{signature(x = "ContactGroup")}: Method to
      transform each contact matrix to logCPM scale.}
  }
}
\section{Utilities}{
  \describe{
    \item{}{\code{contacts(x)}, \code{contacts(x) <- value}: Method to
      extract the list of contact matrices from a
      ContactGroup. \code{value} is expected to be a list object.}
    \item{distanceIdx}{\code{signature(cg = "ContactGroup",
        threshold="ANY", step="ANY")}: Method to identify which matrix bands
      are no more than \code{threshold} bins apart, where each bin
      represents \code{step} base pairs.}
  }
}

\value{
  A ContactGroup object.
}

\details{
  The ContactGroup class contains a set of contact matrices in the slot
  `contacts`.  All matrices are required to be of the same
  dimensionality.  `ContactGroup()` expects a list of symmetric matrices to
  be passed to the contstructor.  Data about these contact matrices is
  held in two other slots.  Data about the genomic loci represented in
  the ContactGroup is found in the `rowData` slot as a GenomicRanges objects, and sample-level
  information is located in the `colData` slot as a DataFrame.
}

\references{
  Law, C.W., Chen, Y., Shi, W. and Smyth G.K.
  \emph{voom: Precision weights unlock linear model analysis tools for RNA-seq read counts.}
  Genome Biology 2014, 15:R29. \doi{10.1186/gb-2014-15-2-r29}.
}


\examples{
data(cgEx)

cgEx[1,]
cgEx[,1]

cd <- colData(cgEx)
colData(cgEx) <- cd

gr <- rowData(cgEx)
rowData(cgEx) <- gr

cgEx

cl <- contacts(cgEx)
contacts(cgEx) <- cl

d.idx <- distanceIdx(cgEx, 1e7, 4e4)

libs <- librarySize(cgEx)

cgEx.cpm <- logCPM(cgEx)

## below, upper.mats.list is a list of upper triangular matrices
## SampleData is a DataFrame of sample data and LociData is a GenomicRanges objec
\dontrun{
  MatsList <- lapply(upper.mats.list, function(M) M[lower.tri(M)] = M[upper.tri(M)])
  cg <- ContactGroup(LociData, MatsList, SampleData)
}

}
\keyword{classes}

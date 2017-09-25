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
\section{Constructor}{
  \describe{
    \item{}{
      \code{ContactGroup(rowData = GRanges, contacts = list, colData = DataFrame)}:
      Creates a ContactGroup object.
      \describe{
        \item{\code{rowData}:}{Object of class \code{GenomicRanges} equal
          in length to the number of rows/columns in contact matrices. }
        \item{\code{contacts}:}{Object of class \code{list} that contains
          all contact matrices.}
        \item{\code{colData}:}{Object of class \code{DataFrame} containing
          sample-level information.}
      }
    }
  }
}

\section{Methods}{
  In the code snippets below, \code{x} is a ContactGroup object.
  \describe{
    \item{[}{\code{signature(x = "ContactGroup", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{colData}{\code{signature(x = "ContactGroup")}: ... }
    \item{colData<-}{\code{signature(x = "ContactGroup", value = "DataFrame")}: ... }
    \item{dim}{\code{signature(x = "ContactGroup")}: ... }
    \item{rowData}{\code{signature(x = "ContactGroup")}: ... }
    \item{rowData<-}{\code{signature(x = "ContactGroup")}: ... }
    \item{show}{\code{signature(object = "ContactGroup")}: ... }
  }
}
\section{Utilities}{
  \describe{
    \item{[}{\code{signature(x = "ContactGroup", i = "ANY", j = "ANY",
	drop = "ANY")}: ... }
    \item{}{\code{x[i = "ANY", j = "ANY",
        drop = "ANY"]}:}{
      Allows for subsetting \code{x} through use
      of \code{i} or of samples through \code{j}.
    }
    \item{}{\code{colData(x)}, \code{colData(x) <- value}:}{
      Get or set sample-level information about samples in 
      \code{x}.  \code{value} is expected to be a \code{DataFrame} object.
    }
    \item{}{\code{dim(x)}:}{
      Obtain the dimensions of a ContactGroup.  Returns 2 values: one
      representing the number of loci and another representing the number of samples.
    }
    \item{}{\code{rowData(x)}, \code{rowData(x) <- value}:}{
      Get or set the GenomicRanges object describing the
      loci in the ContactGroup.  \code{value} is expected to be a
      GenomicRanges object.
    }
    \item{}{\code{show(x)}:}{
      Method to display summary information about a ContactGroup: the
      number of bins, the width of the bins and the number of samples.
    }
    \item{}{contacts(x)}, \code{contacts(x) <- value}:}{
    Method to extract the list of contact matrices from a
    ContactGroup. \code{value} is expected to be a list object.
  }
  }
}
\references{
%%  ~~put references to the literature/web site here~~
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

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

\title{Class \code{"ContactGroup"}}
\description{
  The ContactGroup class represents a collection of contact matrices
  which represent observations from different samples on the same set of
  genomic loci.
}
\usage{
\S4method{ContactGroup}{GRanges, list, DataFrame}{rowData, contacts, colData}
}

\section{Slots}{
  \describe{
    \item{\code{rowData}:}{Object of class \code{"GenomicRanges"} equal
  in length to the number of rows/columns in contact matrices. }
    \item{\code{contacts}:}{Object of class \code{"list"} that contains
  all contact matrices.}
    \item{\code{colData}:}{Object of class \code{"DataFrame"} containing
  sample-level information.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "ContactGroup", i = "ANY", j = "ANY",
        drop = "ANY")}:}{
      Allows for subsetting contact matrices through use
      of \code{i} or of samples through \code{j}.
    }
    \item{colData}{\code{signature(x = "ContactGroup")}:}{
      Accessor that retrieves sample-level information about samples in ContactGroup.
    }
    \item{colData<-}{\code{signature(x = "ContactGroup", value =
        "DataFrame")}:}{
      Set the \code{colData} of a ContactGroup.
    }
    \item{dim}{\code{signature(x = "ContactGroup")}:}{
      Obtain the dimensions of a ContactGroup.  Returns 1 value
  representing the number of loci and another representing the number of samples.
    }
    \item{rowData}{\code{signature(x = "ContactGroup")}:}{
      Accessor that retrieves the GenomicRanges object describing the
  loci in the ContactGroup.
    }
    \item{rowData<-}{\code{signature(x = "ContactGroup")}:}{
      Set the \code{rowData} of a ContactGroup.
    }
    \item{show}{\code{signature(object = "ContactGroup")}:}{
      Method to display summary information about a ContactGroup: the
  number of bins, the width of the bins and the number of samples.
    }
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

}
\keyword{classes}

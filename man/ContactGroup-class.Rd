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
  A class for holding Hi-C contact matrices for multiple samples.
}
\usage{
## Constructor
ContactGroup()

## Accessors
colData()
rowData()
contact()

## Utilities
}
\arguments{
  \item{}{}
}

\section{Slots}{
  \describe{
    \item{\code{rowData}:}{Object of class \code{"GenomicRanges"} ~~ }
    \item{\code{contacts}:}{Object of class \code{"list"} ~~ }
    \item{\code{colData}:}{Object of class \code{"DataFrame"} ~~ }
  }
}
\section{Methods}{
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
  }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\examples{
showClass("ContactGroup")
}
\keyword{classes}

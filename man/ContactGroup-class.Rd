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
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ContactGroup", ...)}.
%%  ~~ describe objects here ~~ 
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
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("ContactGroup")
}
\keyword{classes}

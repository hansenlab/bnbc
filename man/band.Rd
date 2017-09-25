\name{band}
\alias{band}
\title{
  Get Band  
}
\description{
  Get or set band from matrix.
}

\section{Methods}{
  \describe{
    \item{band}{\code{signature(x = "matrix", band.no = "ANY")}:}{
      Extracts band from matrix.
    }
    \item{band<-}{\code{signature(x = "matrix", band.no = "ANY", value =
        "ANY")}:}{
      Set the value of a matrix band to \code{value}.
    }
  }
}

\arguments{
  \item{x}{A matrix.}
  \item{band.no}{Integer specifying which matrix band.  \code{band.no} = 1 retrieves
    the main diagonal.}
  \item{value}{A scalar or vector equal in length to the matrix band.}
    
}

\value{
  A matrix band.
}

\details{
  A matrix band is the set of elements in a matrix from a specific off-diagonal.
}
\examples{
data(chr22)
tact.1 <- contacts(chr22)[[1]]
b2 <- band(tact.1, 2)
band(tact.1, 2) <- b2
}
\seealso{
  \link{ContactGroup}
  \link{getBandIdx}
  \link{band}
  
}

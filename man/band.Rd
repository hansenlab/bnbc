\name{band}
\alias{band}
\alias{band<-}
\title{
  Get Band  
}
\description{
  Get or set band from matrix.
}

\usage{
band(x, band.no)
band(x, band.no) <- value
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
data(cgEx)
tact.1 <- contacts(cgEx)[[1]]
b2 <- band(tact.1, 2)
band(tact.1, 2) <- b2
}
\seealso{
  \link{ContactGroup},
  \link{getBandIdx},
  \link{band}
  
}

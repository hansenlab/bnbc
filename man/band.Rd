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
band(mat, band.no)
band(mat, band.no) <- value
}

\arguments{
  \item{mat}{A matrix.}
  \item{band.no}{Integer specifying which matrix band.  \code{band.no} = 1 retrieves
    the main diagonal.}
  \item{value}{A scalar or vector equal in length to the matrix band.}
    
}

\value{
  A matrix band in the form of a vector.
}

\details{
  A matrix band is the set of elements in a matrix from a specific off-diagonal.
}

\seealso{
  \link{getBandIdx}
}

\examples{
mat <- matrix(1:9, 3, 3)
band(mat, band.no = 2)
mat
band(mat, band.no = 2) <- c(9,10)
mat

data(cgEx)
tact.1 <- contacts(cgEx)[[1]]
b2 <- band(tact.1, 2)
band(tact.1, 2) <- b2
}

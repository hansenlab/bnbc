\name{getBandMatrix}
\alias{getBandMatrix}
\title{
  Get Band Matrix 
}
\description{
  Get band matrix from ContactGroup.
}
\usage{getBandMatrix(cg, band.no=1)}

\arguments{
  \item{cg}{A ContactGroup object.}
  \item{band.no}{Integer specifying which matrix band.  \code{band.no} = 1 retrieves
  the main diagonal.}
}

\value{
  A matrix with one column per sample in the ContactGroup and number of
  rows equal to the length of the matrix band.
}

\details{
  A band matrix is a matrix whose columns are the \code{band.no}-th
  off-diagonal of each sample's contact matrix.  If there are \code{k}
  samples and matrix band \code{band.no} has \code{r} entries, then the
  returned band matrix is of dimension \code{r} x \code{k}.
}
\examples{
data(cgEx)
b2 <- getBandMatrix(cgEx, 2)
}
\seealso{
  \link{ContactGroup},
  \link{getBandIdx},
  \link{band}
}

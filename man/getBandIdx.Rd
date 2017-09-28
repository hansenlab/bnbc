\name{getBandIdx}
\alias{getBandIdx}
\title{
  Get Band Indices
}
\description{
  Get the indices corresponding to a matrix band.
}
\usage{getBandIdx(n, band.no)}

\arguments{
  \item{n}{The number of rows/columns of a contact matrix}
  \item{band.no}{Integer specifying which matrix band.  \code{band.no} = 1 retrieves
  the main diagonal.}
}

\value{
  A matrix with 2 columns and as many rows as entries in the matrix band.
}

\details{
  This function is used in subsetting contact matrices, primarily in
  \link{getBandMatrix}.  However, users wishing to extract band matrices
  directly may find this useful
}
\examples{
data(cgEx)
b2.idx <- getBandIdx(nrow(cgEx), 2)
}
\seealso{
  \link{ContactGroup},
  \link{getBandMatrix},
  \link{band}
  
}

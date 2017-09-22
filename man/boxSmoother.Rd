\name{boxSmoother}
\alias{boxSmoother}
\docType{data}
\title{
  Using Box Smoother
}
\description{
  This function applies a square smoothing kernel of radius \code{h} to
  all contact matrices in a ContactGroup object.
}
\usage{boxSmoother(cg, h, mc.cores)}

\arguments{
  \item{cg}{A ContactGroup object.}
  \item{h}{The desired smoother radius.  This is an integer.}
  \item{mc.cores}{The number of cores to be used.}
}

\value{
  A ContactGroup object is returned that contains the smoothed matrices.
}

\details{
  This method applies a square smoothing kernel of radius \code{h} to
  all contact matrices in a ContactGroup object.  Specifying radius
  \code{h} implies that the width of the kernel is \eqn{2*h + 1} matrix cells.
}
\examples{
data(chr22)
chr22.smooth <- boxSmoother(chr22, h=5, mc.core=1)
}
\seealso{
  \link{ContactGroup}
}
\keyword{datasets}

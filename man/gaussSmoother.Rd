\name{gaussSmoother}
\alias{gaussSmoother}
\title{
  Using Gaussian Smoother
}
\description{
  This function applies a square Gaussain smoothing kernel of width
  \code{radius} with standard deviation \code{sigma} to
  all contact matrices in a ContactGroup object.
}
\usage{boxSmoother(cg, radius, sigma)}

\arguments{
  \item{cg}{A ContactGroup object.}
  \item{radius}{The desired smoother radius.  This is an integer.}
  \item{sigma}{The standard deviation of the smoother.}
}

\value{
  A ContactGroup object is returned that contains the smoothed matrices.
}

\details{
  This method applies a square smoothing kernel of width \code{radius}
  with standard deviation \code{sigma} to all contact matrices in a ContactGroup object.  
}
\examples{
data(chr22)
chr22.smooth <- gaussSmoother(chr22, radius=3, sigma=0.5)
}
\seealso{
  \link{ContactGroup}
}


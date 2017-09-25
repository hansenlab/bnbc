\name{smoothing}
\alias{boxSmoother}
\alias{gaussSmoother}
\title{
  Smoothing Operations
}
\description{
  These functions apply a smoothing kernel to all contact matrices in a
  ContactGroup object.
}
\section{Methods}{
  In the code snippets below, \code{cg} is a ContactGroup object.
  \describe{
    \item{}{
      \code{boxSmoother(cg, h, mc.cores)}:
      Apply a box smoother.
    }
    \item{}{
      \code{gaussSmoother(cg, radius, sigma, mc.cores)}:
      Apply a Gaussian smoother}
  }
}
\arguments{
  \item{cg}{A ContactGroup object.}
  \item{h}{The desired smoother radius.  Only applies to box smoother. This is an integer.}
  \item{radius}{The desired smoother width. Only applies to Gaussian
    smoother.  This is an integer.}
  \item{sigma}{The desired smoother standard deviation. Only applies to
    Gaussian smoother.  This is a positive number.}
  \item{mc.cores}{The number of cores to be used.}
}

\value{
  A ContactGroup object is returned that contains the smoothed matrices.
}

\details{
  \code{boxSmoother} applies a square smoothing kernel of radius \code{h} to
  all contact matrices in a ContactGroup object.  Specifying radius
  \code{h} implies that the width of the kernel is \eqn{2*h + 1} matrix
  cells.

  \code{gaussSmoother} applies a square Gaussain smoothing kernel of width
  \code{radius} with standard deviation \code{sigma} to
  all contact matrices in a ContactGroup object.
}
\examples{
data(cgEx)
cgEx.smooth <- boxSmoother(cgEx, h=5, mc.cores=1)
cgEx.smooth <- gaussSmoother(cgEx, radius=3, sigma=0.5, mc.cores=1)
}
\seealso{
  \link{ContactGroup}
}


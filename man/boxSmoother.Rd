\name{boxSmoother}
\alias{boxSmoother}
\alias{gaussSmoother}
\title{
  2D smoothers for Hi-C contact matrices.
}
\description{
  2D smoothers for Hi-C contact matrices.
}
\usage{
boxSmoother(cg, h, mc.cores=1)
gaussSmoother(cg, radius=3, sigma=0.5)
}
\arguments{
  \item{cg}{}
  \item{h}{}
  \item{mc.cores}{}
  \item{radius}{}
  \item{mc.sigma}{}
}
\details{}
\value{
  An object of class \code{ContactGroup}.
}
\references{
}
\examples{
data(chr22)
chr22 <- boxSmoother(chr22, h = 5)
}

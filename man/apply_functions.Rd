\name{cgApply}
\alias{cgApply}
\alias{cgBandApply}
\title{
  Apply-type methods
}
\description{
  These functions are apply-type functions for ContactGroup objects.
}
\usage{
cgApply(cg, FUN, mc.cores=1, ...)
cgBandApply(cg, FUN, nbands=NULL, mc.cores=1, bstart=2, ...)
}
\arguments{
  \item{cg}{A ContactGroup object.}
  \item{FUN}{A function to be applied. For \code{cgApply} this function
    should operate on a square matrix. For \code{cgBandApply} this
    function should operate on a band, ie. a vector (see \code{\link{band}}).}
  \item{mc.cores}{The number of cores to be used.  Defaults to 1}
  \item{bstart}{The first band to apply a function to.  Defaults to 2.
    Only applicable to \code{cgBandApply}.}
  \item{nbands}{The last band to apply a function to.  Default is
    \code{nrow(cg) - 1}. Only applicable to \code{cgBandApply}.}
  \item{...}{Passed to \code{mclapply}.}
}
\value{
  For \code{cgApply}, a ContactGroup object.  For \code{cgBandApply}, a
  list whose elements are the returned value of \code{FUN}.
}
\details{
  These methods make it easy to apply functions to either all contact
  matrices or a set of bands in all contact matrices.  Both methods
  accept a function \code{FUN}.  For \code{cgApply}, the first argument
  should be \code{cg}, the contact group itself.  For
  \code{cgBandApply}, the first argument should also be \code{cg}, and
  the second argument should be a specific band number.  Additionally,
  the bands to be iterated are specified through \code{bstart:nbands}: \code{bstart}
  indicates the starting band, and \code{nbands} indicates the
  last band.
}
\seealso{
  \code{\link{ContactGroup}},
  \code{\link{getBandMatrix}},
  \code{\link{band}}
}
\examples{
data(cgEx)
cgEx.1 <- cgApply(cgEx, FUN=function(xx){ xx + 1 })
band.matrix.list <- cgBandApply(cgEx, FUN=getBandMatrix, bstart=2, nbands=5)
}


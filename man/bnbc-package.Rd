\name{bnbc-package}
\alias{bnbc-package}
\docType{package}
\title{
\packageTitle{bnbc}
}
\description{
\packageDescription{bnbc}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bnbc}
\packageIndices{bnbc}

The package implements the bnbc method for normalizing Hi-C data across
samples. The name is short for band-wise normalization and batch
correction. The main workhorse is the \code{bnbc} function. We recommend
using smoothing and library size normalization first.

The package implements the \code{ContactGroup} class for storing
multiple Hi-C contact matrices. This is most naturally done with one
object per chromosome, which is ugly.

We also have functions for applying over a \code{ContactGroup}
(\code{cgApply}) and working with matrix bands \code{band, getBandIdx}.

}
\author{
\packageAuthor{bnbc}

Maintainer: \packageMaintainer{bnbc}
}
\references{
  Fletez-Brant et al.
  \emph{Distance-dependent between-sample normalization for Hi-C
  experiments}. In preparation.
}
\keyword{ package }
\seealso{
  \code{\link{bnbc}}, \code{\link{ContactGroup}}, \code{\link{band}},
  \code{\link{cgApply}}. 
}
\examples{
data(cgEx)
batches <- colData(cgEx)$Batch
cgEx.cpm <- logCPM(cgEx)
cgEx.smooth <- boxSmoother(cgEx, 5, mc.cores=1)
cgEx.bnbc <- bnbc(cgEx.smooth, batches, 1e7, 4e4, bstart=2, nbands=4)
}

\name{bnbc}
\alias{bnbc}
\title{
  Normalize Contact Matrices with BNBC
}
\description{
  Applies BNBC method to normalize contact matrices.
}
\usage{bnbc(cg, batch, threshold = NULL, step = NULL, qn = TRUE, nbands
= NULL, mod = NULL, mean.only = FALSE, tol = 5, bstart = 2)}

\arguments{
  \item{cg}{A ContactGroup object.}
  \item{batch}{A single batch indicator variable.}
  \item{threshold}{The maximum distance interacting loci are allowed to
    be separated by.}
  \item{step}{The step size, or the number of bases a contact matrix
    cell represents.}
  \item{qn}{Whether to apply quantile normalization on each band matrix.
    Defaults to TRUE.}
  \item{bstart}{The first band to normalize. Defaults to 2.}
  \item{nbands}{The last band to normalize. Defaults to \code{nrow(cg) -
      1}.}
  \item{mod}{A model matrix specifying which sample information is to be
    preserved by ComBat.  Optional.}
  \item{mean.only}{Whether ComBat should not correct for batch effect in the
    variances of band matrix rows. Defaults to FALSE, which means
    variances are corrected.  Set to TRUE if there is only one observation
    per batch.}
  \item{tol}{The number of significant digits for which the mean value
    of a band matrix must be greater than 0 to be processed by ComBat.}
}

\value{
  A ContactGroup object for which matrix bands \code{bstart:nbands} have
  had BNBC applied.
}

\details{
  This function applies BNBC to the set of contact matrices and returns
  a ContactGroup object with matrix bands \code{bstart:nbands}
  corrected.  For those rows in the matrix bands  which cannot be
  corrected we set all elements to 0.

}
\references{
  Johnson, W.E., Li, C. and Rabinovic, A. (2007)
  Adjusting batch effects in microarray expression data using empirical
  Bayes methods.
  \emph{Biostatistics}, \bold{8}.
}

\examples{
data(cgEx)
batches <- colData(cgEx)$Batch
cgEx.cpm <- logCPM(cgEx)
cgEx.bnbc <- bnbc(cgEx.cpm, batches, 1e7, 4e4, bstart=2, nbands=4)
}
\seealso{
  \link{ContactGroup}
  \link{logCPM}
}

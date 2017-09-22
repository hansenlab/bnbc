\name{logCPM}
\alias{logCPM}
\title{
  LogCPM Transform Contact Matrices
}
\description{
  Apply \code{logCPM} transformation to each contact matrix in ContactGroup.
}
\usage{logCPM(cg)}

\arguments{
  \item{x}{A ContactGroup object.}
}

\value{
  A ContactGroup whose contact matrices are now on the logCPM scale.
}

\details{
  The \code{logCPM} transformation for a given contact matrix is defined
  to be \eqn{logCPM(x) = \log\bigg(\frac{x + 0.5}{L + 1}\bigg)}, where
  \code{L} is the library size of the contact matrix.
}
\references{
  Law, C.W., Chen, Y., Shi, W. and Smyth G.K. (2014)
  voom: Precision weights unlock linear model analysis tools for RNA-seq read counts.
  \emph{Genome Biol.}, \bold{15}.
}

\examples{
data(chr22)
chr22.cpm <- logCPM(chr22)
}
\seealso{
  \link{ContactGroup}
  \link{librarySize}
}

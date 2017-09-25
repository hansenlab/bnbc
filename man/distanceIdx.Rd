\name{logCPM}
\alias{logCPM}
\title{
  Interactions at set distance(s).
}
\description{
  Get indices of matrix bands which represent interactions at set distance(s).
}
\usage{distanceIdx(cg, threshold, step)}

\arguments{
  \item{cg}{A ContactGroup object.}
  \item{threshold}{The maximum distance interacting loci are allowed to
    be separated by.}
  \item{step}{The step size, or the number of bases a contact matrix
  cell represents.}
}

\value{
  The set of matrix bands which represent interactions
  bewteen loci separated by a distance less than or equal to \code{threshold}.
}

\details{
  This function returns the set of matrix bands which represent interactions
  bewteen loci separated by a distance less than or equal to \code{threshold}.

}

\examples{
data(chr22)
d.idx <- distanceIdx(chr22, 1e7, 4e4)
}
\seealso{
  \link{ContactGroup}
}

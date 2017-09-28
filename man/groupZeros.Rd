\name{groupZeros}
\alias{getGroupZeros}
\alias{dropGroupZeros}
\title{
  Group Zero Operations
}
\description{
  These functions find and remove rows in the set of contact matrices
  for which all elements in the row are 0 in all samples.
}
\usage{
getGroupZeros(cg)
dropGroupZeros(cg, g0s)
}
\arguments{
  \item{cg}{A ContactGroup object.}
  \item{g0s}{A list of elements identified as group zeros.}
}

\value{
  A ContactGroup object with the group zeros removed from all rows in
  all contact matrices.
}

\details{
  Group zeros are those rows for which all elements of the row in all
  samples are 0.  These can impact estimation of features such as A/B
  compartment status and so should be removed for many analyses.
}
\examples{
data(cgEx)
g0s <- getGroupZeros(cgEx)
cgEx <- dropGroupZeros(cgEx, g0s)
}
\seealso{
  \link{ContactGroup}
}


\name{librarySize}
\alias{librarySize}
\title{
  Get Library Sizes from ContactGroup 
}
\description{
  This function computes the library size of each contact matrix in a
  ContactGroup (see \code{Details}).
}
\usage{librarySize(cg)}

\arguments{
  \item{cg}{A ContactGroup object.}
}

\value{
  A vector library sizes, one for each contact matrix in a ContactGroup.
}

\details{
  The library size of a contact matrix is defined to be the sum of
  entries of theupper triangle of that contact matrix, excluding that matrix's
  diagonal.  This method computes the library size for each contact
  matrix in a ContactGroup.
}
\examples{
data(chr22)
librarySize(chr22)
}
\seealso{
  \link{ContactGroup}
  \link{logCPM}
}

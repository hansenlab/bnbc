\name{cooler_methods}
\alias{getGenomeIdx}
\alias{getChrCGFromCools}
\alias{cg2bedgraph2}

\title{
  Methods for manipulating cooler files
}

\description{
  These are a set of methods for working with data in cooler file
  format. At present this assumes there is one sample's data at one
  resolution per cooler file, for all observations from one sample.
}

\usage{
getGenomeIdx(file, step)
getChrCGFromCools(bin.ixns, files, chr, colData)
cg2bedgraph2(cg, out.dir)
}

\arguments{
  \item{file}{The name of a cooler file.}
  \item{step}{The resolution of the data inside the cooler file.}
  \item{bin.ixns}{A list that contains a list of interactions and the
    bins involved; output of \code{getGenomeIdx}.}
  \item{files}{A vector of cooler file names.}
  \item{chr}{The target chromosome to be read.}
  \item{colData}{A \code{data.frame} or \code{DataFrame} of metadata for the
    \code{ContactGroup} object.}
  \item{cg}{A \code{ContactGroup} object.}
  \item{out.dir}{A directory in which individual bedgraph2 (BG2) files
    are to be written.}
}

\value{
  For \code{getGenomeIdx}, a list that contains a set of all
  interactions, as well as a listing of all bins present in the cooler
  file \code{cool.fh}. For \code{getChrCGFromCools}, a
  \code{ContactGroup} object. There is nothing returned by
  \code{cg2bedgraph2}.
}

\details{
  These methods allow for the normalization of cooler files. Users must
  first generate an index of the genome represented in all cooler
  files (it is assumed the same loci pair's interactions are observed in
  all samples), which allows for subsequent efficient I/O
  operations. \code{getGenomeIdx} generates this index. The index can
  then be used to pull data for one chromosome for all cooler files
  using \code{getChrCGFromCools}. Users can then follow the standard
  pipeline, and save their data in bedgraph2 (BG2) format using
  \code{cg2bedgraph2}. cooler provides a tool to convert this format to cooler
  and users are encouraged to make use of this tool.
}

\examples{
coolerDir <- system.file("cooler", package = "bnbc")
cools <- list.files(coolerDir, pattern="cool$", full.names=TRUE)

step <- 4e4

bin.ixns.list <- bnbc:::getGenomeIdx(cools[1], step)
bin.ixns <- bin.ixns.list$bin.ixns

data(cgEx)
cool.cg <- bnbc:::getChrCGFromCools(bin.ixns,
                                    files = cools,
                                    chr = "chr22",
                                    colData = colData(cgEx)[1:2,])
all.equal(contacts(cgEx)[[1]], contacts(cool.cg)[[1]])
}

\seealso{
  \link{ContactGroup}
}


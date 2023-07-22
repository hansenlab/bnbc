\name{cooler_methods}
\alias{getGenomeIdx}
\alias{getChrCGFromCools}
\alias{cg2bedgraph2}

\title{
  Methods for manipulating cooler files
}

\description{
  These are a set of methods for working with data in cooler file
  format. 
}

\usage{
getChrIdx(chr.length, step)
getChrCGFromCools(files, chr, step, index.gr, work.dir, exp.name, coldata)
cg2bedgraph2(cg, out.dir, prefix)
}

\arguments{
  \item{chr.length}{The length of a chromosome.}
  \item{step}{The resolution of the data inside the cooler file.}
  \item{files}{A vector of cooler file names.}
  \item{chr}{The target chromosome to be read.}
  \item{index.gr}{A \code{GRanges} object, can be output from
    \code{getChrIdx}.}
  \item{work.dir}{Directory for saving temporary files.}
  \item{exp.name}{The name of the experiment, will be appended all
  output file names.}
  \item{coldata}{A \code{data.frame} or \code{DataFrame} of metadata for the
    \code{ContactGroup} object.}
  \item{cg}{A \code{ContactGroup} object.}
  \item{out.dir}{A directory in which individual bedgraph2 (BG2) files
    are to be written.}
  \item{prefix}{A prefix for all output files; e.g. "treatment_study_".}
}

\value{
  For \code{getChrIdx} a \code{GRanges} object with coordinates
  for each bin. For \code{getChrCGFromCools}, a
  \code{ContactGroup} object. There is nothing returned by
  \code{cg2bedgraph2}.
}

\details{
  These methods allow for the normalization of cooler files. Users must
  create their own index, for which we provide \code{getChrIdx}, which
  is an input into \code{getChrCGFromCools}, which uses HiCBricks to
  access the cooler files, and returns a \code{ContactGroup}
  object. Users can then follow the standard pipeline, and save their
  data in bedgraph2 (BG2) format using \code{cg2bedgraph2}. cooler
  provides a tool to convert this format to cooler and users are
  encouraged to make use of this tool. Note that HiCBricks expects
  multiple resolutions in the cooler file.
}

\dontrun{
coolerDir <- system.file("cooler", package = "bnbc")
cools <- list.files(coolerDir, pattern="cool$", full.names=TRUE)

step <- 4e4

ixns <- bnbc:::getGenomeIdx(seqlengths(BSgenome.Hsapiens.UCSC.hg19)["chr22"], step)

data(cgEx)
cool.cg <- bnbc:::getChrCGFromCools(files = cools,
chr = "chr22",
step=step,
index.gr=ixns,
work.dir="tmp.dir",
exp.name="example_case",
colData = colData(cgEx)[1:2,])
all.equal(contacts(cgEx)[[1]], contacts(cool.cg)[[1]])
}

\seealso{
  \link{ContactGroup}
}


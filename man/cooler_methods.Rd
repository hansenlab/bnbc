\name{cooler_methods}
\alias{getGenomeIdx}
\alias{getChrCGFromCools}
\alias{cg2bg2}

\title{
  Methods for manipulating cooler files
}

\description{
  These are a set of methods for working with data in cooler file
  format. At present this assumes there is one sample's data at one
  resolution per cooler file, for all observations from one sample.
}

\usage{
getGenomeIdx(cool.fh, step)
getChrCGFromCools(bin.ixns, cools, tgt.chr, meta.df)
cg2bg2(cg, out.dir)
}

\arguments{
  \item{cool.fh}{A file path pointing to a cooler file}
  \item{step}{The resolution of the data inside the cooler file}
  \item{bin.ixns}{A list that contains a list of interactions and the
    bins involved; output of \code{getGenomeIdx}}
  \item{cools}{A list of cooler file names}
  \item{tgt.chr}{The target chromosome to be normalized}
  \item{meta.df}{A \code{data.frame} of metadata for the
    \code{ContactGroup} object}
  \item{cg}{A \code{ContactGroup} object}
  \item{out.dir}{A directory in which individual BG2 files are to be
    written}
}

\value{
  For \code{getGenomeIdx}, a list that contains a set of all
  interactions, as well as a listing of all bins present in the cooler
  file \code{cool.fh}. For \code{getChrCGFromCools}, a
  \code{ContactGroup} object. There is nothing returned by
  \code{cg2bg2}.
}

\details{
  These methods allow for the normalization of cooler files. Users must
  first generate an index of the genome represented in all cooler
  files (it is assumed the same loci pair's interactions are observed in
  all samples), which allows for subsequent efficient I/O
  operations. \code{getGenomeIdx} generates this index. The index can
  then be used to pull data for one chromosome for all cooler files
  using \code{getChrCGFromCools}. Users can then follow the standard
  pipeline, and save their data in BG2 format using
  \code{cg2bg2}. cooler provides a tool to convert this format to cooler
  and users are encouraged to make use of this tool.
}


\dontrun{
data(cgEx)
step <- 4e4

cools <- list.files(dir, pattern="cool", full.names=TRUE)

bin.ixns.list <- getGenomeIdx(cools[1], step)
bin.ixns <- bin.ixns.list$bin.ixns
bins <- bin.ixns.list$bins

batch <- colData(cgEx)$Batch

cool.cg <- getChrCGFromCools(bin.ixns, cools, "chr22", colData(cgEx))
cool.cg.cpm <- logCPM(cool.cg)
cg.smooth <- boxSmoother(cool.cg.cpm, h=5)

cool.bnbc <- bnbc(cg.smooth, batch=colData(cgEx.smooth)$Batch,
                  threshold=1e7, step=4e4, nbands=11, verbose=FALSE)

}

\seealso{
  \link{ContactGroup}
}
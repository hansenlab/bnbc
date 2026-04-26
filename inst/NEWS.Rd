\name{bnbcnews}
\title{bnbc News}
\encoding{UTF-8}

\section{Version 1.25}{
  \itemize{
    \item Replaced HiCBricks with direct rhdf5 calls in
    \code{getChrCGFromCools} to fix a segfault in HDF5 cleanup that
    was crashing the vignette build. HiCBricks is no longer a
    dependency.
    \item Fixed an issue with build error because of missing Suggests
    \item Added support for cooler files using HiCBricks (this
    happened in an earlier verion).
    \item Various documentation fixes.
  }
}

\section{Version 1.1.2}{
  \itemize{
    \item Removed unneeded so files.
  }
}

\section{Version 0.99}{
  \itemize{
    \item Initial release to Bioconductor.
  }
}

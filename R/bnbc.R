library(SummarizedExperiment)
library(GenomicRanges)
library(Matrix)
library(preprocessCore)
library(sva)

hasRcpp <- TRUE

## do roxygen documentation eventually

tryCatch({
  library(Rcpp)
}, error=function(e){
  message("Rcpp not found, cannot use C-based update\n")
  hasRcpp <- FALSE
})

personSub <- function(se, persons){
    idx <- which(metadata(se)$people == persons)
    se2 <- se
    assay(se2, "tacts") <- assay(se2, "tacts")[,,idx]
    se2
}

getBandIdx <- function(n, band.no){
    ## number of elements is number of rows - (band.no - 1)
    n.elems <- n - (band.no - 1)
    ## number of rows is the number of elements
    i <- 1:n.elems
    ## columns are the last n.elem columns, starting at band.no
    j <- band.no:n
    ## matrix allows for direct usage as indices
    return(cbind(i=i, j=j))
}

idxSwap <- function(bandidx){
    nbandidx <- bandidx
    nbandidx[,1] <- bandidx[,2]
    nbandidx[,2] <- bandidx[,1]
    nbandidx
}

getBandMatrix <- function(cg, band.no=1){
  arr <- contacts(cg)
  nrows <- nrow(cg)
  nppl <- length(arr)
  tmp <- do.call(cbind, lapply(1:nppl, function(ii){
    arr[[ii]][getBandIdx(nrows, band.no)] }))
  colnames(tmp) <- rownames(metadata(cg))
  tmp
}

if(hasRcpp){
  sourceCpp("src/update_bands.cpp") ## much faster

  bnbcC <- function(cg, threshold, step, batch,
                    wts.cg=NULL, qn=TRUE, nbands=NULL, mod=NULL,
                    mean.only=FALSE, tol=5, bstart=2){
    dims <- dim(cg)
    nppl <- dims[3]
    tacts <- contacts(cg)
    if (!is.null(wts.cg)){
      wts.tacts <- contacts(wts.cg)
    }
    if(is.null(nbands)){
      nbands <- distanceIdx(cg, threshold, step)
    }
    wts.mat <- NULL
    mat.list <- list()
    for (ii in bstart:nbands){
      if (ii %% 50 == 0){ cat(".") }
      mat <- getBandMatrix(tacts, nrow(tacts[[1]]), ii, nppl)
      mat.good <- 1:nrow(mat)
      if (qn){ mat <- normalize.quantiles(mat, FALSE) }
      if(!mean.only){
        batchvars <- bandLevelBatchVars(mat, batch)
        mat.good <- abs(rowMeans(mat)) > 0 & round(rowMeans(batchvars), tol)  > 0
      }
      if(!is.null(wts.cg)){
        wts.mat <- getBandMatrix(wts.tacts, nrow(tacts[[1]]), ii, nppl)
      }
      ## mat[mat.good,] <- shutUpComBat(mat[mat.good,], wts=wts.mat, batch)
      ## mat[!mat.good,] <- 0
      matFix <- function(mat, mat.good, batch , wts.mat, mod, mean.only){
        mat <- tryCatch({
          ## mat[mat.0s,] <- shutUpComBat(mat[mat.0s,], batch)
          mat[mat.good,] <- shutUpComBat(mat[mat.good,], wts=wts.mat, batch, mod=mod,
                                         mean.only=mean.only)
          mat[!mat.good,] <- 0
          mat
        }, error=function(e){
          print(ii)
          mat <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
          mat
        })
      } 
      mat <- matFix(mat, mat.good, batch, wts.mat, mod, mean.only)
      tacts <- updateBand(tact_list=tacts, idx=getBandIdx(nrow(tacts[[1]]), ii)-1,
                          band=mat)
    }
    make.sym <- function(mat){
      mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
      mat
    }
    tacts <- lapply(tacts, make.sym)
    new.cg <- cg
    contacts(new.cg) <- tacts
    new.cg
  }
}

bnbc <- function(cg, threshold, step, batch,
                 wts.cg=NULL, qn=TRUE, nbands=NULL, mod=NULL,
                 mean.only=FALSE, tol=5, bstart=2){
    dims <- dim(cg)
    nppl <- dims[3]
    tacts <- contacts(cg)
    if (!is.null(wts.cg)){
      wts.tacts <- contacts(wts.cg)
    }
    if(is.null(nbands)){
        nbands <- distanceIdx(cg, threshold, step)
    }
    wts.mat <- NULL
    for (ii in bstart:nbands){
      if (ii %% 50 == 0){ cat(".") }
      mat <- getBandMatrix(tacts, nrow(tacts[[1]]), ii, nppl)
      mat.good <- 1:nrow(mat)
      if (qn){ mat <- normalize.quantiles(mat, FALSE) }
      if(!mean.only){
        batchvars <- bandLevelBatchVars(mat, batch)
        mat.good <- abs(rowMeans(mat)) > 0 & round(rowMeans(batchvars), tol)  > 0
      }
      if(!is.null(wts.cg)){
        wts.mat <- getBandMatrix(wts.tacts, nrow(tacts[[1]]), ii, nppl)
      }
      ## mat[mat.good,] <- shutUpComBat(mat[mat.good,], wts=wts.mat, batch)
      ## mat[!mat.good,] <- 0
      matFix <- function(mat, mat.good, batch , wts.mat, mod, mean.only){
        mat <- tryCatch({
          ## mat[mat.0s,] <- shutUpComBat(mat[mat.0s,], batch)
          mat[mat.good,] <- shutUpComBat(mat[mat.good,], wts=wts.mat, batch, mod=mod,
                                       mean.only=mean.only)
          mat[!mat.good,] <- 0
          mat
        }, error=function(e){
          print(ii)
          mat <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
          mat
        })
        mat
      } ## ugh 
      mat <- matFix(mat, mat.good, batch, wts.mat, mod, mean.only)
      for (jj in 1:ncol(mat)){
            band(tacts[[jj]], ii) <- mat[,jj]
      }
    }
    contacts(cg) <- tacts
    cg
}

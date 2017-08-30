library(SummarizedExperiment)
library(GenomicRanges)
library(Matrix)
library(preprocessCore)
library(sva)
library(EBImage)

## do roxygen documentation eventually

smoothCG <- function(cg, radius=3, sigma=0.5) capply(cg, gblur,
                                                     sigma=sigma, radius=radius,
                                                     boundary="replicate")

personSub <- function(se, persons){
    idx <- which(metadata(se)$people == persons)
    se2 <- se
    assay(se2, "tacts") <- assay(se2, "tacts")[,,idx]
    se2
}

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

smoothCG <- function(cg, radius=3, sigma=0.5) {
    capply(cg, gblur, sigma=sigma, radius=radius, boundary="replicate")
}

personSub <- function(se, persons){
    idx <- which(metadata(se)$people == persons)
    se2 <- se
    assay(se2, "tacts") <- assay(se2, "tacts")[,,idx]
    se2
}

make.sym <- function(mat){
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    mat
}

bnbc <- function(cg, batch, threshold=NULL, step=NULL,
                  qn=TRUE, nbands=NULL, mod=NULL,
                  mean.only=FALSE, tol=5, bstart=2){
    dims <- dim(cg)
    nppl <- dims[3]
    tacts <- contacts(cg)
    if(is.null(nbands)){
        nbands <- distanceIdx(cg, threshold, step)
    }
    mat.list <- list()
    for (ii in bstart:nbands){
        if (ii %% 50 == 0){ cat(".") }
        mat <- getBandMatrix(cg, ii)
        mat.good <- 1:nrow(mat)
        if (qn){ mat <- normalize.quantiles(mat, FALSE) }
        if(!mean.only){
            batchvars <- bandLevelBatchVars(mat, batch)
            mat.good <- abs(rowMeans(mat)) > 0 & round(rowMeans(batchvars), tol)  > 0
        }
        mat[mat.good,] <- ComBat(mat[mat.good,], batch, mod=mod,
                                 mean.only=mean.only)
        mat[!mat.good,] <- 0
        tacts <- updateBand(tact_list=tacts, idx=getBandIdx(nrow(tacts[[1]]), ii)-1,
                            band=mat)
    }
    tacts <- lapply(tacts, make.sym)
    new.cg <- cg
    contacts(new.cg) <- tacts
    new.cg
}

## bnbc <- function(cg, batch, threshold=NULL, step=NULL,
##                  qn=TRUE, nbands=NULL, mod=NULL,
##                  mean.only=FALSE, tol=5, bstart=2){
##     dims <- dim(cg)
##     nppl <- dims[3]
##     tacts <- contacts(cg)
##     if(is.null(nbands)){
##         nbands <- distanceIdx(cg, threshold, step)
##     }
##     for (ii in bstart:nbands){
##         if (ii %% 50 == 0){ cat(".") }
##         mat <- getBandMatrix(cg, ii)
##         mat.good <- 1:nrow(mat)
##         if (qn){ mat <- normalize.quantiles(mat, FALSE) }
##         if(!mean.only){
##             batchvars <- bandLevelBatchVars(mat, batch)
##             mat.good <- abs(rowMeans(mat)) > 0 & round(rowMeans(batchvars), tol)  > 0
##         }
##         mat[mat.good,] <- ComBat(mat[mat.good,], batch, mod=mod,
##                                  mean.only=mean.only)
##         mat[!mat.good,] <- 0
##         for (jj in 1:ncol(mat)){
##             band(tacts[[jj]], ii) <- mat[,jj]
##         }
##     }
##     contacts(cg) <- tacts
##     cg
## }

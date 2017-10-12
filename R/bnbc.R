bnbc <- function(cg, batch, threshold=NULL, step=NULL,
                  qn=TRUE, nbands=NULL, mod=NULL,
                  mean.only=FALSE, tol=5, bstart=2){
    nppl <- ncol(cg)
    tacts <- contacts(cg)
    if(is.null(nbands)){
        nbands <- distanceIdx(cg, threshold = threshold, step = step)
    }
    mat.list <- list()
    for (ii in bstart:nbands){
        if (ii %% 50 == 0){ cat(".") }
        mat <- getBandMatrix(cg, ii)
        mat.good <- seq_len(nrow(mat))
        if (qn) {
            mat <- normalize.quantiles(mat, copy = FALSE)
        }
        if(!mean.only){
            batchvars <- bandLevelBatchVars(mat, batch)
            mat.good <- abs(rowMeans(mat)) > 0 &
                round(rowMeans(batchvars), tol)  > 0
        }
        tryCatch({
            suppressMessages({
                mat[mat.good,] <- ComBat(mat[mat.good,], batch, mod=mod,
                                         mean.only=mean.only)
            })
        }, error=function(e){ warning() })
        mat[!mat.good,] <- 0
        tacts <- updateBand(tact_list=tacts,
                            idx=getBandIdx(nrow(tacts[[1]]), ii)-1,
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

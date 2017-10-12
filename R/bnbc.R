bnbc <- function(cg, batch, threshold=NULL, step=NULL,
                  qn=TRUE, nbands=NULL, mod=NULL,
                  mean.only=FALSE, tol=5, bstart=2, verbose = TRUE){
    nppl <- ncol(cg)
    tacts <- contacts(cg)
    if(is.null(nbands)){
        nbands <- distanceIdx(cg, threshold = threshold, step = step)
    }
    stopifnot(bstart >= 1, bstart <= nbands)
    mat.list <- list()
    if(verbose) pb <- txtProgressBar(style = 3)
    for (ii in bstart:nbands){
        if (verbose && ii %% 50 == 0){ setTxtProgressBar(pb, (ii - bstart) / nbands) }
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
        }, error=warning)
        mat[!mat.good,] <- 0
        tacts <- updateBand(tact_list=tacts,
                            idx=getBandIdx(nrow(tacts[[1]]), ii)-1,
                            band=mat)
    }
    if(verbose) close(pb)
    make.sym <- function(mat){
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        mat
    }
    tacts <- lapply(tacts, make.sym)
    new.cg <- cg
    contacts(new.cg) <- tacts
    new.cg
}

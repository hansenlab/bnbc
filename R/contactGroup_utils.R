logCPM <- function(x){
    stopifnot(is(x, "ContactGroup"))
    libs <- (librarySize(x) + 1) / 10^6
    contacts(x) <- lapply(seq_len(ncol(x)), function(ii) {
        log(contacts(x)[[ii]] + 0.5 / libs[ii] + 1)
    })
    x
}

librarySize <- function(x) {
    stopifnot(is(x, "ContactGroup"))
    sapply(contacts(x), function(mat) {
        sum(mat[upper.tri(mat)])
    })
}

gaussSmoother <- function(cg, radius=3, sigma=0.5) {
    stopifnot(is(cg, "ContactGroup"))
    cgApply(cg, gblur, sigma=sigma, radius=radius, boundary="replicate")
}

boxSmoother <- function(cg, h, mc.cores=1){
    stopifnot(is(cg, "ContactGroup"))
    smoother <- function(xx, h){
        size <- (2 * h + 1)
        sbox <- makeBrush(size=size, shape="box")/size^2
        filter2(xx, sbox, boundary="replicate")
    }
    cgApply(cg, smoother, h=h, mc.cores=mc.cores)
}

cgBandApply <- function(cg, FUN, nbands=NULL, mc.cores=1, bstart=2, ...){
    stopifnot(is(cg, "ContactGroup"))
    if(is.null(nbands)){
        nbands <- nrow(cg) - 1
    }
    if (mc.cores == 1){
        return(sapply(bstart:nbands, function(ii){
            FUN(cg, ii, ...)
        }))
    }
    if(mc.cores > 1){
        return(mclapply(bstart:nbands, function(ii){
            FUN(cg, ii, ...)
        }, mc.cores=mc.cores))
    }
}

cgApply <- function(cg, FUN, mc.cores=1, ...){
    stopifnot(is(cg, "ContactGroup"))
    if (mc.cores == 1) {
        contacts(cg) <- lapply(contacts(cg), FUN = FUN, ...)
    }
    if(mc.cores > 1) {
        contacts(cg) <- mclapply(contacts(cg), FUN, ..., mc.cores = mc.cores)
    }
    cg
}

distanceIdx <- function(cg, threshold, step){
    stopifnot(is(cg, "ContactGroup"))
    return((1:nrow(cg))[which(1:nrow(cg) * step <= threshold)])
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
    is(cg, "ContactGroup")
    arr <- contacts(cg)
    idxs <- getBandIdx(nrow(cg), band.no = band.no)
    tmp <- do.call(cbind, lapply(seq_len(ncol(cg)), function(ii){
                              arr[[ii]][idxs]
                          }))
    colnames(tmp) <- rownames(colData(cg))
    tmp
}

bandLevelBatchVars <- function(mat, batches){
    sapply(unique(batches), function(ii){
        rowVars(mat[,batches == ii])
    })
}

band <- function(mat, band.no){
    return(mat[getBandIdx(nrow(mat), band.no)])
}

`band<-` <- function(mat, band.no, value){
    b.idx <- getBandIdx(nrow(mat), band.no)
    mat[b.idx] <- value
    mat[idxSwap(b.idx)] <- value
    mat
}

## the set of contact matrices have rows/columns of g0s removed
## as they are only 0.  also remove corresponding ranges in granges ('gr').
## first add a new toplevel list element which is a granges object
## it has colData col specifying whether a range was kept/no
## this list element has no accessor method, so it does not disrupt
## any existing code (cf validContactGroup, which only checks for specific items)
## use should occur after !is.null(cg$g0s) test
dropGroupZeros <- function(cg, g0s){
    stopifnot(is(cg, "ContactGroup"))
    cg$g0s <- rowData(cg)
    mcols(cg$g0s)$kept <- TRUE
    mcols(cg$g0s)$kept[g0s] <- FALSE
    cg <- cg[-g0s]
    cg
}

getGroupZeros <- function(cg){
    stopifnot(is(cg, "ContactGroup"))
    Reduce(intersect,
           lapply(contacts(cg), function(xx){
               xx[getBandIdx(nrow(xx), 1)] <- 0
               return(which(rowMeans(xx) == 0))
           })
           )
}
utils::globalVariables(c("chrom", "chrom.1", "bin1_id", "bin2_id"))

logCPM <- function(x){
    stopifnot(is(x, "ContactGroup"))
    libs <- (librarySize(x) + 1) / 10^6
    contacts(x) <- lapply(seq_len(ncol(x)), function(ii) {
        log((contacts(x)[[ii]] + 0.5) / libs[ii])
    })
    x
}

librarySize <- function(x) {
    stopifnot(is(x, "ContactGroup"))
    sapply(contacts(x), function(mat) {
        sum(mat[upper.tri(mat)])
    })
}

gaussSmoother <- function(cg, radius=3, sigma=0.5, mc.cores=1) {
    stopifnot(is(cg, "ContactGroup"))
    cgApply(cg, gblur, mc.cores=mc.cores, sigma=sigma,
            radius=radius, boundary="replicate")
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
    stopifnot(bstart >= 1, bstart <= nbands)
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
    return(seq_len(nrow(cg))[which(seq_len(nrow(cg)) * step <= threshold)])
}

getBandIdx <- function(n, band.no){
    stopifnot(band.no >= 1, n >= band.no)
    ## number of elements is number of rows - (band.no - 1)
    n.elems <- n - (band.no - 1)
    ## number of rows is the number of elements
    i <- seq_len(n.elems)
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

getGenomeIdx <- function(file, step){
    cool.h <- H5Fopen(file)
    bins <- data.frame(cool.h$bins)
    ixns <- data.frame(cool.h$pixels[1:2])
    H5Fclose(cool.h)
    bin.ixns <- data.frame(ixns,
                       bins[ixns[, 1]+1,],
                       bins[ixns[, 2]+1,])
    rownames(bin.ixns) <- rownames(ixns)
    for (ii in c(1, 2, 4, 5, 7,8)){
        bin.ixns[[ii]] <- as.numeric(bin.ixns[[ii]])
    }
    bin.ixns$dist <- NA
    cis.dists <- bin.ixns[bin.ixns$chrom == bin.ixns$chrom.1, "start.1"] -
        bin.ixns[bin.ixns$chrom == bin.ixns$chrom.1, "start"]
    bin.ixns[bin.ixns$chrom == bin.ixns$chrom.1,]$dist <- cis.dists
    bin.ixns <- data.table(bin.ixns, keep.rownames=TRUE)
    bin.ixns$i <- bin.ixns$start/step
    bin.ixns$j <- bin.ixns$start.1/step
    setkey(bin.ixns, bin1_id, bin2_id)
    list(bins=bins, bin.ixns=bin.ixns)
}

getChrCGFromCools <- function(bin.ixns, files, chr, colData){
    bin.ixns.now <- bin.ixns[chrom == chr & chrom.1 == chr]
    i1 <- bin.ixns.now$i + 1
    j1 <- bin.ixns.now$j + 1
    tgt.rows <- as.numeric(bin.ixns.now$rn)
    max.ij <- max(max(i1), max(j1))

    cool.dat <- lapply(files, function(xx){
        mat.df <- h5read(xx, name="pixels/count", index=list(tgt.rows))
        mat <- matrix(0, nrow=max.ij, ncol=max.ij)
        mat[cbind(i1, j1)] <- mat.df
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]        
        mat
    }) 
    gr <- GRanges(unique(unique(bin.ixns.now[, c("chrom", "end", "start")])))
    ContactGroup(rowData=gr, contacts=cool.dat, colData=colData)
}

cg2bedgraph2 <- function(cg, out.dir){
    chr.ref <- data.frame(rowData(cg))
    ijs <- which(upper.tri(contacts(cg)[[1]], diag=TRUE), arr.ind=TRUE)
    devnull <- lapply(seq_along(contacts(cg)), function(ii){
        fname <- file.path(out.dir, paste0(paste(data.frame(colData(cg)[ii,]),
                                                 collapse="_"), ".bg2"))
        message(fname)
        out.df <- data.frame(chr.ref[ijs[,1], 1:3],
                             chr.ref[ijs[,2], 1:3],
                             contacts(cg)[[ii]][ijs])
        colnames(out.df) <- c("chrom1", "start1", "end1", 
                              "chrom2", "start2", "end2",
                              "count")
        write.table(out.df, file=fname, row.names=FALSE, col.names=FALSE, sep="\t",
                    quote=FALSE)
    })
    message("completed")
}


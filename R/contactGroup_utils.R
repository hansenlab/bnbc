librarySize <- function(x){
  libs <- sapply(1:dim(x)[3], function(ii){
    tmp <- contacts(x)[[ii]]
    sum(tmp[upper.tri(tmp)])
  })
}

setMethod("librarySize", "ContactGroup",
          function(x){
            libs <- sapply(1:dim(x)[3], function(ii){
              tmp <- contacts(x)[[ii]]
              sum(tmp[upper.tri(tmp)])
            })
            libs
          })


librarySizeScale <- function(x, k=1e6, a=0.5, b=1){
    libs <- librarySize(x)
    contacts(x) <- lapply(1:dim(x)[3], function(ii){
        return(log((contacts(x)[[ii]] + a )/(libs[ii] + b) * k))
    })
    x
}


bapply <- function(cg, FUN, nbands=NULL, ncores=1, bstart=2, ...){
  if(is.null(nbands)){
    nbands <- nrow(cg) - 1
  }
  if (ncores == 1){
    return(sapply(bstart:nbands, function(ii){
      FUN(cg, ii, ...)
    }))
  }
  else{
    return(mclapply(bstart:nbands, function(ii){
      FUN(cg, ii, ...)
    }, mc.cores=ncores))
  }
}

capply <- function(cg, FUN, ncores=1, ...){
  cg2 <- cg
  if (ncores == 1){
    contacts(cg2) <- lapply(contacts(cg), function(xx){
      FUN(xx, ...)
    })
  }
  else{
    contacts(cg2) <- mclapply(contacts(cg), function(xx){
      FUN(xx, ...)
    }, mc.cores=ncores)
  }
  cg2
}

distanceIdx <- function(cg, threshold, step){
    return(c(1:dim(cg)[1])[which(1:dim(cg)[1] * step <= threshold)])
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
  colnames(tmp) <- rownames(colData(cg))
  tmp
}

bandLevelBatchVars <- function(mat, batches){
  sapply(unique(batches), function(ii){
    rowVars(mat[,batches == ii])
  })
}

band <- function(mat, band.no){ return(mat[getBandIdx(nrow(mat), band.no)]) }

assign("band<-",
       function(mat, band.no, value){
           b.idx <- getBandIdx(nrow(mat), band.no)
           mat[b.idx] <- value
           mat[idxSwap(b.idx)] <- value
           mat
           })

## the set of contact matrices have rows/columns of g0s removed
## as they are only 0.  also remove corresponding ranges in granges ('gr').
## first add a new toplevel list element which is a granges object
## it has colData col specifying whether a range was kept/no
## this list element has no accessor method, so it does not disrupt
## any existing code (cf validContactGroup, which only checks for specific items)
## use should occur after !is.null(cg$g0s) test
dropGroupZeros <- function(cg, g0s){
  cg$g0s <- granges(cg)
  mcols(cg$g0s)$kept <- TRUE
  mcols(cg$g0s)$kept[g0s] <- FALSE
  cg <- cg[-g0s]
  cg
}

getGroupZeros <- function(cg){
  Reduce(intersect, lapply(contacts(cg), function(xx){
    xx[getBandIdx(nrow(xx), 1)] <- 0
    return(which(rowMeans(xx) == 0)) }))
}

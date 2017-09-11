setOldClass("ContactGroup")

ContactGroup <- function(gr, contacts, metadata){
    out <- list(gr = gr, contacts = contacts, metadata = metadata)
    class(out) <- "ContactGroup"
    stopifnot(validContactGroup(out))
    out
}

## by convention, sample dim should be 3rd array dimension
validContactGroup <- function(x) {
    stopifnot(is.list(x))
    stopifnot(class(x) == "ContactGroup")
    stopifnot(all(c("gr", "contacts") %in% names(x)))
    stopifnot(is(x$gr, "GRanges"))
    stopifnot(length(x$gr) == dim(x$contacts[[1]]))
    nrowcol <- t(sapply(x$contacts, function(xx){
        return(list(nrow(xx), ncol(xx)))
    }))
    stopifnot(all(nrowcol[1,1] == unlist(nrowcol[1,])))
    stopifnot(all(nrowcol[1,1] == unlist(nrowcol[,1])))
    TRUE
}

print.ContactGroup <- function(x, ...) {
    cat("Genome contact matrix with\n")
    cat(sprintf(" %s bins\n", length(x$gr)))
    cat(sprintf(" %s kb average width per bin\n", round(mean(width(x$gr)) / 1000)))
    cat(sprintf(" %s samples\n", length(x$contacts)))
    invisible(x)
}


metadata <- function(x){
    x$metadata
}

assign("metadata<-",
       function(x, value){
           x$metadata <- value
           x
       })

setMethod("granges", "ContactGroup",
          function(x, use.mcols = FALSE, ...){
              x$gr
          })

setMethod("mcols", "ContactGroup",
          function(x, use.mcols = FALSE, ...){
              mcols(x$gr)
          })


contacts <- function(x){
    x$contacts
}

assign("contacts<-",
       function(x, value){
           x$contacts <- value
           x
       })

assign("granges<-",
       function(x, value){
           x$gr <- value
           x
       })


length.ContactGroup <- function(x){
    length(x$gr)
}

dropPeople <- function(x, idx){
    contacts(x) <- contacts(x)[idx]
    metadata(x) <- metadata(x)[idx,]
    x
}

## gets i bins
assign("[.ContactGroup",
       function(object, i){
           object$gr <- object$gr[i]
           object$contacts <- lapply(object$contacts, function(xx){
             return(xx[i, i, drop = FALSE]) })
             ## return(xx[i, drop = FALSE]) })
           object
       })

setMethod("start", "ContactGroup",
          function(x, ...) {
              start(granges(x), ...)
          })

setMethod("end", "ContactGroup",
          function(x, ...) {
              end(granges(x), ...)
          })

setMethod("width", "ContactGroup",
          function(x) {
              width(granges(x))
          })

setMethod("seqnames", "ContactGroup",
          function(x) {
              seqnames(granges(x))
          })

setMethod("seqlevels", "ContactGroup",
          function(x) {
              seqlevels(granges(x))
          })

setMethod("seqlengths", "ContactGroup",
          function(x) {
              seqlengths(granges(x))
          })

setMethod("seqinfo", "ContactGroup",
          function(x) {
              seqinfo(granges(x))
          })

setMethod("genome", "ContactGroup",
          function(x) {
              genome(granges(x))
          })

librarySize <- function(x){
  libs <- sapply(1:dim(x)[3], function(ii){
    tmp <- contacts(x)[[ii]]
    sum(triu(tmp))
  })
}

setMethod("librarySize", "ContactGroup",
          function(x){
            libs <- sapply(1:dim(x)[3], function(ii){
              tmp <- contacts(x)[[ii]]
              sum(triu(tmp))
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


dim.ContactGroup <- function(x){
    c(dim(contacts(x)[[1]]), length(contacts(x)))
}

nrow.ContactGroup <- function(x){
    dim(contacts(x)[[1]])[1]
}

ncol.ContactGroup <- function(x){
    dim(contacts(x)[[1]])[2]
}

setMethod("subsetByOverlaps",
signature(query = "ContactGroup", subject = "GenomicRanges"),
function(query, subject, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end", "within", "equal"),
         ignore.strand = FALSE, ...) {
    ov <- findOverlaps(query = granges(query), subject = subject,
                       maxgap = maxgap, minoverlap = minoverlap,
                       type = match.arg(type), select = "first",
                       ignore.strand = ignore.strand, ... )
    query[!is.na(ov)]
})

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "GenomicRanges"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = granges(query), subject = subject,
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "CollapsedVCF"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = granges(query), subject = granges(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "CollapsedVCF", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = granges(query), subject = granges(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = granges(query), subject = granges(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "GenomicRanges", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = query, subject = granges(subject),
                           maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
})

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
  colnames(tmp) <- rownames(metadata(cg))
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
## it has metadata col specifying whether a range was kept/no
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

setClass("ContactGroup",
         rowData = "GenomicRanges",
         contacts = "list",
         colData = "DataFrame")


setValidity("ContactGroup", function(object) {
    txt <- NULL
    contactTxt <- "all components of `contacts` has to be square matrices of the same dimension"
    if(!all(sapply(x@contacts, is.matrix)))
        txt <- contactTxt
    nrowcol <- sapply(x@contacts, dim)
    if(all(nrowcol[1,1] == nrowcol[,1]))
        txt <- contactTxt
    if(all(nrowcol[1,] == nrowcol[2,]))
        txt <- contactTxt
    if(nrow(x@colData) != length(x@contacts))
        txt <- c(txt, "the length of `contacts` has to be the same as the number of rows of `colData`")
    if(!all(names(x@contacts) == rownames(x@colData)))
        txt <- c(txt, "the names of `contacts` has to be equal to the rownames of `colData`" <- )
    if(length(x@rowData) != nrow(x@contacts[[1]]))
        txt <- c(txt, "the length of `rowData` should be equal to the number of rows of the matrices in `contacts`")
    txt
})

ContactGroup <- function(rowData, contacts, colData){
    out <- new("ContactGroup", rowData = rowData, contacts = contacts, colData = colData)
    out
}

## validContactGroup <- function(x) {
##     ## by convention, sample dim should be 3rd array dimension
##     stopifnot(class(x) == "ContactGroup") ##
##     stopifnot(is.list(x))##
##     stopifnot(all(c("rowData", "contacts") %in% names(x)))##
##     stopifnot(is(x$rowData, "GRanges"))##
##     stopifnot(length(x$rowData) == dim(x$contacts[[1]]))
##     nrowcol <- t(sapply(x$contacts, function(xx){
##         return(list(nrow(xx), ncol(xx)))
##     }))
##     stopifnot(all(nrowcol[1,1] == unlist(nrowcol[1,])))
##     stopifnot(all(nrowcol[1,1] == unlist(nrowcol[,1])))
##     TRUE
## }

print.ContactGroup <- function(x, ...) {
    cat("Genome contact matrix with\n")
    cat(sprintf(" %s bins\n", length(x$rowData)))
    cat(sprintf(" %s kb average width per bin\n", round(mean(width(x$rowData)) / 1000)))
    cat(sprintf(" %s samples\n", length(x$contacts)))
    invisible(x)
}

setMethod("colData", signature(x = "ContactGroup"),
          function(x, ...) x$colData)

setReplaceMethod("colData", signature(x = "ContactGroup", value = "ANY"),
                 function(x, ..., value) {
    x$colData <- value
})

setMethod("rowData", signature(x = "ContactGroup"),
          function(x, use.names = TRUE, use.mcols = TRUE) {
    x$rowData
})

setReplaceMethod("rowData", signature(x = "ContactGroup", value = "ANY"),
                 function(x, value){
    x$rowData <- value
    x
})

setMethod("mcols", "ContactGroup",
          function(x, use.mcols = FALSE, ...){
              mcols(rowData(x))
          })

## setMethod("start", "ContactGroup",
##           function(x, ...) {
##               start(granges(x), ...)
##           })

## setMethod("end", "ContactGroup",
##           function(x, ...) {
##               end(granges(x), ...)
##           })

## setMethod("width", "ContactGroup",
##           function(x) {
##               width(granges(x))
##           })

## setMethod("seqnames", "ContactGroup",
##           function(x) {
##               seqnames(granges(x))
##           })

## setMethod("seqlevels", "ContactGroup",
##           function(x) {
##               seqlevels(granges(x))
##           })

## setMethod("seqlengths", "ContactGroup",
##           function(x) {
##               seqlengths(granges(x))
##           })

## setMethod("seqinfo", "ContactGroup",
##           function(x) {
##               seqinfo(granges(x))
##           })

## setMethod("genome", "ContactGroup",
##           function(x) {
##               genome(granges(x))
##           })

setMethod("ncol", signature(x = "ContactGroup"),
          function(x) {
    ncol(contacts(x)[[1]])
})

setMethod("nrow", signature(x = "ContactGroup"),
          function(x) {
    nrow(contacts(x)[[1]])
})

setMethod("dim", signature(x = "ContactGroup"),
          function(x) {
    c(dim(contacts(x)[[1]]), length(contacts(x)))
})


contacts <- function(x){
    x$contacts
}

assign("contacts<-",
       function(x, value){
           x$contacts <- value
           x
       })

dropPeople <- function(x, idx){
    contacts(x) <- contacts(x)[idx]
    colData(x) <- colData(x)[idx,]
    x
}

## gets i bins
assign("[.ContactGroup",
       function(object, i, j){
           object$rowData <- object$rowData[i]
           object$contacts <- lapply(object$contacts[j], function(xx){
               return(xx[i, i, drop = FALSE]) })
           object$colData <- object$colData[j,]
           object
})

## setMethod("subsetByOverlaps",
## signature(query = "ContactGroup", subject = "GenomicRanges"),
## function(query, subject, maxgap = 0L, minoverlap = 1L,
##          type = c("any", "start", "end", "within", "equal"),
##          ignore.strand = FALSE, ...) {
##     ov <- findOverlaps(query = granges(query), subject = subject,
##                        maxgap = maxgap, minoverlap = minoverlap,
##                        type = match.arg(type), select = "first",
##                        ignore.strand = ignore.strand, ... )
##     query[!is.na(ov)]
## })

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "GenomicRanges"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = rowData(query), subject = subject,
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "CollapsedVCF"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = rowData(query), subject = granges(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "CollapsedVCF", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = granges(query), subject = rowData(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "ContactGroup", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = rowData(query), subject = rowData(subject),
                 maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                 ignore.strand = ignore.strand, ...)
})

setMethod("findOverlaps",
signature(query = "GenomicRanges", subject = "ContactGroup"),
function (query, subject, maxgap = 0L, minoverlap = 1L,
          type = c("any", "start", "end", "within", "equal"),
          select = c("all", "first"), ignore.strand = FALSE, ...) {
    findOverlaps(query = query, subject = rowData(subject),
                           maxgap = maxgap, minoverlap = minoverlap,
                 type = match.arg(type), select = match.arg(select),
                           ignore.strand = ignore.strand, ...)
})


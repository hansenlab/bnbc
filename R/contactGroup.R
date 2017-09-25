setClass("ContactGroup",
         representation(rowData = "GenomicRanges",
                        contacts = "list",
                        colData = "DataFrame"))


setValidity("ContactGroup", function(object) {
    txt <- NULL
    contactTxt <- "all components of `contacts` has to be square matrices of the same dimension"
    if(!all(sapply(object@contacts, is.matrix)))
        txt <- contactTxt
    nrowcol <- sapply(object@contacts, dim)
    if(any(nrowcol[1,1] != nrowcol[,1]))
        txt <- contactTxt
    if(any(nrowcol[1,] != nrowcol[2,]))
        txt <- contactTxt
    if(nrow(object@colData) != length(object@contacts))
        txt <- c(txt, "the length of `contacts` has to be the same as the number of rows of `colData`")
    if(!all(names(object@contacts) == rownames(object@colData)))
        txt <- c(txt, "the names of `contacts` has to be equal to the rownames of `colData`")
    if(length(object@rowData) != nrow(object@contacts[[1]]))
        txt <- c(txt, "the length of `rowData` should be equal to the number of rows of the matrices in `contacts`")
    txt
})

ContactGroup <- function(rowData, contacts, colData){
    ## FIXME: check that the constructor returns a valid object when called as ContactGroup()
    out <- new("ContactGroup", rowData = rowData, contacts = contacts, colData = colData)
    out
}

setMethod("show", signature(object = "ContactGroup"),
          function(object) {
    cat("Object of class `ContactGroup` representing contact matrices with\n")
    cat(sprintf(" %s bins\n", length(rowData(object))))
    cat(sprintf(" %s kb average width per bin\n", round(mean(width(rowData(object))) / 1000)))
    cat(sprintf(" %s samples\n", nrow(colData(object))))
    invisible(object)
})

setMethod("colData", signature(x = "ContactGroup"),
          function(x, ...) {
    x@colData
})

setReplaceMethod("colData", signature(x = "ContactGroup", value = "DataFrame"),
                 function(x, ..., value) {
    x@colData <- value
    x
})

setMethod("rowData", signature(x = "ContactGroup"),
          function(x, ...) {
    x@rowData
})

setReplaceMethod("rowData", signature(x = "ContactGroup"),
                 function(x, ..., value){
    x@rowData <- value
    x
})

contacts <- function(x){
    x@contacts
}

`contacts<-` <- function(x, value){
    x@contacts <- value
    x
}

setMethod("dim", signature(x = "ContactGroup"),
          function(x) {
    c(length(x@rowData), nrow(x@colData))
})

setMethod("[", signature(x = "ContactGroup", i = "ANY", j = "ANY"),
          function(x, i, j, ...) {
    if (!missing(i)) {
        x@rowData <- x@rowData[i]
        x@contacts <- lapply(x@contacts[j], function(xx){
            return(xx[i, i, drop = FALSE]) })
    }
    if(!missing(j)) {
        x@colData <- x@colData[j,]
    }
    x
})

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


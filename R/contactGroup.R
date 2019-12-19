setClass("ContactGroup",
         representation(rowData = "GenomicRanges",
                        contacts = "list",
                        colData = "DataFrame"))


setValidity("ContactGroup", function(object) {
    txt <- NULL
    contactTxt <- "all components of `contacts` has to be square matrices of the same dimension"
    if(!all(sapply(contacts(object), is.matrix)))
        txt <- contactTxt
    nrowcol <- sapply(contacts(object), dim)
    ## special case of CG with no entries
    ## nrowcol is empty list instead of 2d object
    if (is.list(nrowcol)){
      if (length(rowData(object)) == 0 & nrow(colData(object)) == 0){
        return(txt)
      }
    }
    if(any(nrowcol[1,1] != nrowcol[,1]))
        txt <- contactTxt
    if(any(nrowcol[1,] != nrowcol[2,]))
        txt <- contactTxt
    if(nrow(colData(object)) != length(contacts(object)))
        txt <- c(txt, "the length of `contacts` has to be the same as the number of rows of `colData`")
    if(!all(names(contacts(object)) == rownames(colData(object))))
        txt <- c(txt, "the names of `contacts` has to be equal to the rownames of `colData`")
    if(length(rowData(object)) != nrow(contacts(object)[[1]]))
        txt <- c(txt, "the length of `rowData` should be equal to the number of rows of the matrices in `contacts`")
    txt
})

ContactGroup <- function(rowData=GRanges(), contacts=list(), colData=DataFrame()){
    out <- new("ContactGroup", rowData = rowData, contacts = contacts,
               colData = colData)
    out
}

setMethod("show", signature(object = "ContactGroup"),
          function(object) {
    cat("Object of class `ContactGroup` representing contact matrices with\n")
    cat(sprintf(" %s bins\n", length(rowData(object))))
    xx <- round(mean(width(rowData(object))) / 1000)
    cat(sprintf(" %s kb average width per bin\n",
                ifelse(is.na(xx), 0, xx)))
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
    c(length(rowData(x)), nrow(colData(x)))
})

setMethod("[", signature(x = "ContactGroup", i = "ANY", j = "ANY"),
          function(x, i, j, ...) {
    if (!missing(i)) {
        rowData(x) <- rowData(x)[i]
        contacts(x) <- lapply(contacts(x)[j], function(xx){
            return(xx[i, i, drop = FALSE]) })
    }
    if(!missing(j)) {
        colData(x) <- colData(x)[j,]
        contacts(x) <- contacts(x)[j]
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

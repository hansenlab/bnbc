# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

replace <- function(mat, idx, update) {
    .Call(`_bnbc_replace`, mat, idx, update)
}

updateBand <- function(tact_list, idx, band) {
    .Call(`_bnbc_updateBand`, tact_list, idx, band)
}

getBandIdxC <- function(n, band_no) {
    .Call(`_bnbc_getBandIdxC`, n, band_no)
}

getSeq <- function(a, b) {
    .Call(`_bnbc_getSeq`, a, b)
}

updateBands <- function(tact_list, n, bstart, nbands, band_mats) {
    .Call(`_bnbc_updateBands`, tact_list, n, bstart, nbands, band_mats)
}


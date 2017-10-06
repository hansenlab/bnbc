test_bands <- function() {
    mat <- matrix(1:9, 3, 3)
    checkEquals(band(mat, band.no = 2), c(4,8))
    band(mat, band.no = 2) <- c(9,10)
    checkEquals(mat[1,2], 9)
    checkEquals(mat[2,1], 9)
    checkEquals(mat[2,3], 10)
    checkEquals(mat[3,2], 10)
}


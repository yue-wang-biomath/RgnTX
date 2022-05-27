#' getFormatCorrect
#' @description This function makes sure the two input region sets are in the correct format required by RgnTX evaluation functions.
#'
#' @usage getFormatCorrect(A, B)
#'
#' @param A Region set 1. A Granges or GRangesList object.
#' @param B Region set 2. A Granges or GRangesList object.
#'
#' @return A \code{list} object.

getFormatCorrect <- function(A, B){
    if (is(A, "GRangesList")) {
        A <- GRangesList2GRanges(A)
    }
    if (is(B, "GRangesList")) {
        B <- GRangesList2GRanges(B)
    }

    if (!is(A, "GRanges")) {
        stop("A should be either GRanges or GRangesList.")
    }
    if (!is(B, "GRanges")) {
        stop("B should be either GRanges or GRangesList.")
    }

    if (length(A$group) == 0) {
        A$group <- seq_len(length(A))
    }

    if (length(B$group) == 0) {
        B$group <- seq_len(length(B))
    }
    return(list(A,B))
}

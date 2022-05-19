#' getFormatCorrect
#' @export getFormatCorrect
#' @description This function makes sure the two input region sets are in the correct format required by RgnTX evaluation functions.
#'
#' @usage getFormatCorrect(A, B)
#'
#' @param A Region set 1. A Granges or GRangesList object.
#' @param B Region set 2. A Granges or GRangesList object.
#'
#' @return A \code{list} object.
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' A <- randomizeTx(
#'     txdb, trans.ids,
#'     random_num = 20,
#'     random_length = 100
#' )
#' B <- randomizeTx(
#'     txdb, trans.ids,
#'     random_num = 20,
#'     random_length = 100
#' )
#' getFormatCorrect(A, B)

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

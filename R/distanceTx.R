#' Evaluation function
#' @export distanceTx
#' @importFrom IRanges distanceToNearest shift
#' @importFrom methods is
#'
#' @description Evaluation function. This function calculates the mean of the distance from each region of set RS1 to the closest region in RS2.
#'
#' @usage distanceTx(A, B, beta = 0.2, ...)
#'
#' @param A Region set 1. A Granges or GRangesList object.
#' @param B Region set 2. A Granges or GRangesList object.
#' @param beta It is a user-defined argument that can filter out the corresponding percent of largest distance values. Default value is 0.2.
#' @param ... Any additional parameters needed.
#'
#' @return A \code{numeric} object.
#'
#' @seealso \code{\link{overlapWidthTx}}, \code{\link{overlapCountsTx}}
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
#' distanceTx(A, B, beta = 0.2)
distanceTx <- function(A, B, beta = 0.2, ...) {
    A_B <- getFormatCorrect(A, B)
    A <- A_B[[1]]
    B <- A_B[[2]]
    suppressWarnings(
        distance.frame <- data.frame(distanceToNearest(A, B))
    )

    distances <- distance.frame[, "distance"]
    distances <- sort(distances)
    index <- round((1 - beta) * length(distances))
    distances <- distances[seq_len(index)]
    return(mean(distances))
}

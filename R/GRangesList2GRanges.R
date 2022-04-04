#' Convert a GRangesList object to a GRanges object.
#' @export GRangesList2GRanges
#'
#' @description Convert a \code{GRangesList} object to a \code{GRanges} object.
#'
#' @usage GRangesList2GRanges(A = NULL)
#'
#' @param A A \code{GRangesList} object.
#'
#' @return A \code{GRanges} object.
#'
#' @seealso \code{\link{GRanges2GRangesList}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' RS1 <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 100)
#'
#' RS1 <- GRangesList2GRanges(RS1)
GRangesList2GRanges <-
    function(A = NULL) {
        A.df <- data.frame(A)
        if (length(which(is.na(A.df[, "group_name"]) == FALSE)) != 0) {
            R <- GRanges(
                seqnames = A.df[, "seqnames"],
                IRanges(
                    start = A.df[, "start"],
                    end = A.df[, "end"]
                ),
                strand = A.df[, "strand"],
                transcriptsHits = A.df[, "group_name"],
                group = A.df[, "group"]
            )
        } else {
            R <- GRanges(
                seqnames = A.df[, "seqnames"],
                IRanges(
                    start = A.df[, "start"],
                    end = A.df[, "end"]
                ),
                strand = A.df[, "strand"],
                group = A.df[, "group"]
            )
        }
        return(R)
    }

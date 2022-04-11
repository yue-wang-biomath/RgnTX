#' Evaluation function
#' @export overlapCountsTx
#' @importFrom IRanges countOverlaps
#' @importFrom IRanges findOverlaps
#'
#' @description This function receives two region sets and returns the number of their overlaps.
#'
#' @usage overlapCountsTx(A, B, count_once = TRUE, over_trans = TRUE, ...)
#'
#' @param A Region set 1. A \code{GRangesList} object.
#' @param B Region set 2. A \code{GRangesList} object.
#' @param count_once Whether the overlap of multiple B regions with a single A region should be counted once or multiple times.
#' @param over_trans Whether the overlapping is counted over the transcriptome or over the genome.
#' @param ... Any additional parameters needed.
#'
#' @return A \code{numeric} object.
#'
#' @seealso \code{\link{overlapCountsTx}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' exons.tx0 <- exonsBy(txdb)
#' regions.A <- exons.tx0[trans.ids]
#' A <- randomizeTransByOrder(regions.A, random_length = 200)
#' B <- randomizeTransByOrder(regions.A, random_length = 200)
#'
#' overlapCountsTx(A, B)
overlapCountsTx <- function(A, B, count_once = TRUE, over_trans = TRUE, ...) {
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
    suppressWarnings(
        map.df <- data.frame(findOverlaps(A, B))
    )
    if (nrow(map.df) != 0) {
        if (over_trans == TRUE) {
            if (length(A$transcriptsHits) != 0 && length(B$transcriptsHits) != 0) {
                # Filter out overlaps between two regions not on the same transcripts
                A.transcriptsHits <- A$transcriptsHits[map.df[, 1]]
                B.transcriptsHits <- B$transcriptsHits[map.df[, 2]]
                map.df <- map.df[A.transcriptsHits == B.transcriptsHits, ]
            } else {
                print("Either A or B does not provide transcript id information. Their overlapping is counted at genome level rather than transcriptome")
            }
        }

        # Overlapping counts contributed by multiple ranges of the same feature should be only counted once.
        A.group <- A$group[map.df[, 1]]
        B.group <- B$group[map.df[, 2]]
        group.df <- paste0(A.group, "_", B.group)
        group.df.unique <- unique(group.df)

        if (count_once == FALSE) {
            return(length(group.df.unique))
        } else {
            group.num1.end.index <- regexpr("_", group.df.unique)
            df.unique <- substr(group.df.unique, 1, group.num1.end.index - 1)
            return(length(unique(df.unique)))
        }
    } else {
        return(0)
    }
}

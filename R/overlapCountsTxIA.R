#' Evaluation function
#' @export overlapCountsTxIA
#' @importFrom IRanges countOverlaps
#' @importFrom IRanges findOverlaps
#'
#' @description Evaluation function. This function receives a feature set (with isoform ambiguity) and a transcriptome region set (without isoform ambiguity), and returns a weighted number of overlaps between them.
#'
#' @usage overlapCountsTxIA(A, B, ...)
#' @param ... Any additional parameters needed.
#'
#' @param A A feature set, which should be \code{GRanges}.
#' @param B A region set, which should be \code{GRangesList}.
#'
#' @return A \code{numeric} object.
#'
#' @seealso \code{\link{overlapWidthTx}}, \code{\link{distanceTx}}, \code{\link{overlapCountsTx}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' file <- system.file(package="RgnTX", "extdata/m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' RS1 <- m6A_sites_data[1:100]
#'
#' trans.info <- getTransInfo(RS1, txdb)
#' trans.ids <- trans.info[, "trans_ID"]
#'
#' RS2 <- getStopCodon(trans.ids, txdb = txdb)
#'
#' # Evaluation step.
#' orig.ev <- overlapCountsTxIA(RS1, RS2)
overlapCountsTxIA <- function(A, B, ...) {
    trans.info <- getTransInfo(A, txdb)
    trans.id <- trans.info[, "trans_ID"]
    index_features <- trans.info[, "index_features"]
    A_amplify <- A[index_features]
    A_amplify$transcriptsHits <- trans.id
    A_amplify$group <- trans.info[, "index_trans"]
    A_amplify_index <- A_amplify$group
    A_amplify$index_features <- index_features

    A <- A_amplify

    A_transcriptsHits <- as.character(A$transcriptsHits)

    B_names <- names(B)
    intersect_trans <- intersect(A_transcriptsHits, B_names)

    A <- A[is.element(A_transcriptsHits, intersect_trans)]
    A$group <- seq_len(length(A))
    B <- B[is.element(B_names, intersect_trans)]
    A_transcriptsHits <- as.character(A$transcriptsHits)
    B <- B[A_transcriptsHits]

    index_features <- A$index_features
    B <- GRangesList2GRanges(B)
    A$transcriptsHits <- paste0(A$transcriptsHits, "_", A$group)
    B$transcriptsHits <- paste0(B$transcriptsHits, "_", B$group)


    index_weights <- unlist(lapply(index_features, function(x) {
        return(1 / length(which(index_features ==
            x)))
    }))

    suppressWarnings(
    map.df <- data.frame(findOverlaps(A, B))
    )

    A.group <- A$group[map.df[, 1]]
    B.group <- B$group[map.df[, 2]]

    # transcript id information
    if (length(A$transcriptsHits) != 0 & length(B$transcriptsHits) != 0) {
        A.id <- A$transcriptsHits[map.df[, 1]]
        B.id <- B$transcriptsHits[map.df[, 2]]

        map.df <- data.frame(
            A.hit = map.df[, 1], B.hit = map.df[, 2],
            A.id = A.id, B.id = B.id, A.group = A.group, B.group = B.group
        )

        # Filter out wrong mapping Filter out ranges are not on the same transcript.
        map.df <- map.df[A.id == B.id, ]
    }
    A.group <- map.df[, "A.group"]
    B.group <- map.df[, "B.group"]

    # Overlaps over the same group should be only counted once. Filter out repeated counts.
    filterSameGroup <- function(x) {
        A1 <- A.group[x]
        B1 <- B.group[x]
        repeate.index <- which(A.group == A1)
        B.group.rep <- B.group[repeate.index]
        B.group.rep.index <- which(B.group.rep == B1)
        return(repeate.index[min(B.group.rep.index)])
    }

    group.filter <- unlist(lapply(seq_len(nrow(map.df)), filterSameGroup))

    map.df <- map.df[group.filter, ]
    counts.df <- unlist(lapply(map.df[, "A.hit"], function(x) {
        return(index_weights[x])
    }))

    return(sum(counts.df))
}

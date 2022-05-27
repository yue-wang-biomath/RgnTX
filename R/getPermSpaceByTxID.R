#' Get permutation space by specifying transcript ids
#'
#' @export getPermSpaceByTxID
#'
#' @description This function returns 5'UTR/CDS/3'UTR/mRNA/full part of transcriptome regions grouped by corresponding transcript ids.
#'
#' @usage getPermSpaceByTxID(trans_ids = "all", txdb, type = "mature")
#'
#' @param trans_ids A character object. The transcript ids. Default is "all". If it takes the default value "all", the space that users get will be the whole transcriptome.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#'
#' @return A \code{GRangesList} object.
#'
#' @seealso \code{\link{getPermSpaceByType}}, \code{\link{getPermSpaceByFeatures}}
#' @examples
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' permspace <- getPermSpaceByTxID(trans.ids, txdb)
getPermSpaceByTxID <- function(trans_ids = "all", txdb, type = "mature") {

    # make sure the inputs are reasonable
    trans.ids <- trans_ids
    if (!is.character(trans.ids)) {
        stop("trans.ids must be character.")
    }

    regions.A <- getPermSpaceByType(txdb, type = type)

    if (trans.ids[1] == "all") {
        trans.ids <- names(regions.A)
    } else {
        trans.ids <- as.character(trans.ids)
        trans.ids <- trans.ids[is.element(trans.ids, names(regions.A))]
    }

    if (length(trans.ids) != 0) {
        regions.A <- regions.A[trans.ids]
    } else {
        regions.A <- NULL
    }

    return(regions.A)
}

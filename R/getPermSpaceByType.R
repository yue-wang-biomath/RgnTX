#' Get permutation space by specifying type
#' @export getPermSpaceByType
#' @importFrom GenomicFeatures fiveUTRsByTranscript threeUTRsByTranscript cdsBy transcriptsBy exonsBy
#'
#' @description This function can return 5'UTR/CDS/3'UTR/mRNA/full part of transcriptome regions, following the format required by the main permutation test functions.
#'
#' @usage getPermSpaceByType(txdb, type = "mature")
#'
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#'
#' @return A \code{GRangesList} object.
#'
#' @seealso \code{\link{getPermSpaceByTxID}}, \code{\link{getPermSpaceByFeatures}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' permSpace <- getPermSpaceByType(txdb, type = "CDS")
getPermSpaceByType <- function(txdb, type = "mature") {

    # make sure the inputs are reasonable
    if (!is.character(type)) {
        stop("type must be character.")
    }

    regions.A <- NULL
    if (type == "mature") {
        exons.tx0 <- exonsBy(txdb, by = c("tx", "gene"), use.names = FALSE)
        regions.A <- exons.tx0
    }
    if (type == "fiveUTR") {
        fiveUTR.tx0 <- fiveUTRsByTranscript(txdb, use.names = FALSE)
        regions.A <- fiveUTR.tx0
    }
    if (type == "threeUTR") {
        threeUTR.tx0 <- threeUTRsByTranscript(txdb, use.names = FALSE)
        regions.A <- threeUTR.tx0
    }
    if (type == "CDS") {
        cds.tx0 <- cdsBy(txdb, use.names = FALSE)
        regions.A <- cds.tx0
    }
    if (type == "full") {
        trans.tx0 <- transcriptsBy(txdb, "gene")
        trans.tx0.df <- data.frame(trans.tx0)
        trans.tx0.id <- trans.tx0.df[, "tx_id"]
        # RefSeqID
        RefSeqID <- as.character(trans.tx0.id)
        # targetName
        targetName <- trans.tx0.df[, "seqnames"]
        targetName <- as.character(targetName)
        # strand
        strand <- trans.tx0.df[, "strand"]
        strand <- as.character(strand)
        # blockSizes
        blockSizes <- trans.tx0.df[, "width"]

        # targetStart
        targetStart <- trans.tx0.df[, "start"]

        regions.A <- vector2GRangesList(RefSeqID, targetName, strand, blockSizes, targetStart)
    }

    if (is.null(regions.A)) {
        stop("type must come from one of these options: 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR'.")
    }

    return(regions.A)
}

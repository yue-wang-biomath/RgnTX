#' Get permutation space for features
#' @export getPermSpaceByFeatures
#'
#' @description This function returns a default permutation space for features with isoform ambiguity. The default permutation space of a feature is the aggregate of the multiple transcripts it may overlap with. It requires the input feature to be \code{GRanges} format.
#'
#' @usage getPermSpaceByFeatures(features, txdb, type = "mature")
#'
#' @param features A \code{GRanges} object.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of regions over transcriptome.
#'
#' @return
#' A list object, which contains two elements.
#' \itemize{
#' \item \bold{\code{perm.space:}} A \code{GRangesList} object that includes all the transcripts input features may overlap with.
#' \item \bold{\code{index:}} It contains a series of numbers indicating which feature these transcripts are respectively associated with.
#' }
#'
#' @seealso \code{\link{getPermSpaceByTxID}}, \code{\link{getPermSpaceByType}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' file <- system.file(package="RgnTX", "extdata/m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' permSpace <- getPermSpaceByFeatures(features = m6A_sites_data[1:100], txdb)
getPermSpaceByFeatures <- function(features, txdb, type = "mature") {

    # make sure the inputs are reasonable

    if (!is.character(type)) {
        stop("type must be character.")
    }

    regions.A <- getPermSpaceByType(txdb, type = type)
    A <- features

    if (!is(A, "GRanges")) {
        stop("features must be a GRanges object.")
    } else {
        trans_info <- getTransInfo(A, txdb)

        trans_id <- as.character(trans_info[, "trans_ID"])
        index_methyl <- trans_info[, "index_features"]
        index_trans <- trans_info[, "index_trans"]

        regions.A <- regions.A[trans_id]
    }
    return.list <- list(regions.A, index_methyl)
    names(return.list) <- c("perm.space", "index")
    return(return.list)
}

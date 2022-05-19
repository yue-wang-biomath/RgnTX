#' Randomize features into transcriptome
#' @export randomizeFeaturesTx
#'
#' @description Randomize features into transcriptome.
#'
#' @usage randomizeFeaturesTx(RS, txdb, type = "mature", N = 1, ...)
#'
#' @param RS The feature to be randomized. It should be a \code{GRanges} or \code{GRangesList} object.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @param N The number of iterations.
#' @param ... Any additional parameters needed.
#'
#' @return A \code{GRangesList} object. The name of each element is the id of the transcript where the corresponding range is located.
#'
#' @seealso \code{\link{randomizeTransByOrder}}, \code{\link{randomizeFeaturesTxIA}}, \code{\link{randomizeTx}}
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' RS1 <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 100)
#' RS <- randomizeFeaturesTx(RS1, txdb, N = 1)
randomizeFeaturesTx <- function(RS, txdb, type = "mature", N = 1, ...) {
    randomResults.List <- list()

    # three input ways 1. specify a permutation space 2. specify
    # features without IA 3. specify features with IA

    # This function is the second method.

    # If input RS is GRangesList:
    if (is(RS, "CompressedGRangesList")) {
        if (is.character(names(RS))) {
            features <- RS
        }
     }

    # If input RS is GRanges:
    if (is(RS, "GRanges")) {
        metadata <- mcols(RS)
        trans.id <- metadata$transcriptsHits
        if (length(which(is.na(trans.id) == FALSE)) != 0) {
            features <- GRanges2GRangesList(RS)
        } else {
            stop("Transcript id information is missing. It should be provided by its metadata and this metadata column should be named as 'transcriptsHits'.")
        }
    }

    trans.id <- as.character(names(features))
    random.num <- length(features)
    random.length <- sum(width(features))
    regions.A <- getPermSpaceByTxID(trans.id, txdb = txdb, type = type)

    getRandomFeatures <- function(x) {
        randomResults <- randomizeTransByOrder(regions.A, random.length)
        return(randomResults)
    }

    randomResults.List <- lapply(seq_len(N), getRandomFeatures)
    if (N == 1) {return(randomResults.List[[1]])}
    if (N > 1) {return(randomResults.List)}
}

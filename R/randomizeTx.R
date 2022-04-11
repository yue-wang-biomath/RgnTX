#' Get randmized regions over transcriptome
#' @export randomizeTx
#'
#' @description Pick random regions over specified transcripts.
#'
#' @usage randomizeTx(txdb, trans_ids = 'all',
#' random_num = 100, random_length = 20, type = 'mature', N = 1, ...)
#'
#' @param txdb A TxDb object.
#' @param trans_ids The ids of transcripts, which should be a character object. Random regions will be picked from these transcripts. If this argument takes the default value 'all', the scope of picking random regions will be the whole transcriptome.
#' @param random_num The number of regions to be picked.
#' @param random_length The length of regions to be picked.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @param N Randomization times.
#' @param ... Any additional parameters needed.
#'
#' @return A \code{GRangesList} object. The name of each element is the id of the transcript where the corresponding range is located.
#'
#' @seealso \code{\link{randomizeTransByOrder}}, \code{\link{randomizeFeaturesTx}}, \code{\link{randomizeFeaturesTxIA}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' RS1 <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 100)
randomizeTx <- function(txdb, trans_ids = "all", random_num = 100, random_length = 20,
                        type = "mature", N = 1, ...) {
    randomResults.List <- list()

    # three input ways 1. specify a permutation space 2. specify
    # features without IA 3. specify features with IA

    # This function is the first method.

    # make sure the inputs are reasonable

    if (!is.character(trans_ids)) {
        stop("trans_ids must be character.")
    }

    if (!is.numeric(random_num)) {
        stop("random_num must be numeric.")
    }

    if (!is.numeric(random_length)) {
        stop("random_length must be numeric.")
    }

    if (!is.numeric(N)) {
        stop("N must be numeric.")
    }

    perm.space <- getPermSpaceByTxID(trans_ids = trans_ids, txdb, type = type)

    getRandomResultsList <- function(x) {
        trans.id <- sample(names(perm.space), random_num, replace = TRUE)
        regions.A <- perm.space[trans.id]
        randomResults <- randomizeTransByOrder(regions.A, random_length = random_length)
        return(randomResults)
    }
    randomResults.List <- lapply(seq_len(N), getRandomResultsList)
    if (N == 1) {
        return(randomResults.List[[1]])
    }
    if (N > 1) {
        return(randomResults.List)
    }
}

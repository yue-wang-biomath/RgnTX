#' Randomize features into transcriptome.
#' @export randomizeFeaturesTxIA
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicRanges GRanges
#'
#' @description Randomize features into transcriptome, especially for the features that have isoform ambiguity.
#'
#' @usage randomizeFeaturesTxIA(RS, txdb, type = 'mature', N = 1, ...)
#'
#' @param RS The feature being randomized. It should be a GRanges or GRangesList object.
#' @param txdb A txdb object.
#' @param type This argument receives options 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR', with which user can get corresponding types of transcriptome regions.
#' @param N Randomization times.
#' @param ... any additional parameters needed.
#'
#' @return A \code{GRangesList} object. The name of each element is the id of the transcript where the corresponding range is located.
#'
#' @seealso \code{\link{randomizeTransByOrder}}, \code{\link{randomizeFeaturesTx}}, \code{\link{randomizeTx}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' file <- system.file(package="RgnTX", "extdata/m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' RS1 <- m6A_sites_data[1:100]
#' RS <- randomizeFeaturesTxIA(RS1, txdb, N = 1)
randomizeFeaturesTxIA <- function(RS, txdb, type = "mature", N = 1, ...) {
    randomResults.List <- list()

    # three input ways 1. specify a permutation space 2. specify
    # features without IA 3. specify features with IA

    # This function is related to the third method.

    # make sure the inputs are reasonable.
    if (is(RS, "CompressedGRangesList")) {
        RS <- GRangesList2GRanges(RS)
    }

    if (!is(RS, "GRanges")) {
        stop("RS must be GRanges.")

        if (!is.character(type)) {
            stop("type must be character.")
        }

        if (!is.numeric(N)) {
            stop("N must be numeric.")
        }
    }


    features <- RS

    # get permutation space
    perm.list <- getPermSpaceByFeatures(features, txdb, type = type)
    perm.space <- perm.list$perm.space
    index_genomic <- perm.list$index

    # get random.length
    random.length <- width(RS)

    getRandomFeaturesIA <- function(x) {
        trans.index <- lapply(unique(index_genomic), function(x) {
            counts <- which(index_genomic == x)
            if (length(counts) == 1) {
                return(counts)
            } else {
                return(sample(which(index_genomic == x), 1, replace = FALSE))
            }
        })
        trans.index <- unlist(trans.index)
        regions.A <- perm.space[trans.index]
        randomResults <- randomizeTransByOrder(regions.A, random.length)
        return(randomResults)
    }
    randomResults.List <- lapply(seq_len(N), getRandomFeaturesIA)

    if (N == 1) {
        return(randomResults.List[[1]])
    }
    if (N > 1) {
        return(randomResults.List)
    }
}

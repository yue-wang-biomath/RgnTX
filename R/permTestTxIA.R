#' Perform permutation test
#' @export permTestTxIA
#'
#' @description Perform permutation test for evaluating spatial association between some features (with isoform ambiguity) and a region set. It randomizes the features and compares it with the region set to see if there is an association between the features and the region set. The difference between this function and  \code{\link{permTestTx}} is that it is for RNA-related genomic features that have isoform ambiguity, i.e., features that one does not know which transcript they comes from.
#'
#' @usage permTestTxIA(RS1 = NULL,
#'                     RS2 = NULL,
#'                     txdb = NULL,
#'                     type = 'mature',
#'                     ntimes = 50,
#'                     ev_function_1 = overlapCountsTx,
#'                     ev_function_2 = overlapCountsTx,
#'                     pval_z = FALSE,
#'                     ...)
#'
#' @param RS1 The feature set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param RS2 The region set to be compared with. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between \code{RS1} and \code{RS2}. Default is \code{overlapCountsTx}.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in \code{RSL} and \code{RS2}. Default is \code{overlapCountsTx}.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... Any additional parameters needed.
#'
#' @details \code{permTestTxIA} only needs users to input two region sets. It will automatically randomize the first region set into transcriptome.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following items:
#' \itemize{
#' \item \bold{\code{RSL:}} Randomized region sets of \code{RS1}.
#' \item \bold{\code{RS1:}} The feature set to be randomized.
#' \item \bold{\code{RS2:}} The region set to be compared with the feature set.
#' \item \bold{\code{orig.ev:}} The value of overlapping counts between \code{RS1} and \code{RS2}.
#' \item \bold{\code{rand.ev:}} The values of overlapping counts between each element in \code{RSL} and \code{RS2}.
#' \item \bold{\code{pval:}} p-value of the test.
#' \item \bold{\code{zscore:}} Standard score of the test.
#' }
#'
#' @seealso \code{\link{plotPermResults}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' file <- system.file(package="RgnTX", "extdata/m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' RS1 <- m6A_sites_data[1:500]
#' trans.ids <- getTransInfo(RS1, txdb)[, "trans_ID"]
#' RS2 <- getStopCodon(trans.ids, txdb)
#'
#' permTestTx_results <- permTestTxIA(RS1 = RS1, RS2 = RS2,
#'                             txdb = txdb, ntimes = 5)

permTestTxIA <- function(RS1 = NULL, RS2 = NULL, txdb = NULL, type = "mature", ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...) {
    RSL <- randomizeFeaturesTxIA(RS1, txdb, type, N = ntimes)
    permTestTx.results <- permTestTx_customAll(RSL = RSL, RS1 = RS1, RS2 = RS2, txdb = txdb, ev_function_1 = ev_function_1, ev_function_2= ev_function_2, pval_z = pval_z,...)
    return(permTestTx.results)
}

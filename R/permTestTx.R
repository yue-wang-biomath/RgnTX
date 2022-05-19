#' Perform permutation test
#' @export permTestTx
#'
#' @description Perform permutation test for evaluating spatial association between a feature set and a region set.
#'
#' @usage permTestTx(RS1 = NULL, RS2 = NULL, txdb = NULL, type = "mature",
#' ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx,
#' pval_z = FALSE, ...)
#'
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
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
#' exons.tx0 <- exonsBy(txdb)
#' trans.ids <- sample(names(exons.tx0), 500)
#'
#' A <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 100)
#' B <- c(randomizeTx(txdb, trans.ids, random_num = 75, random_length = 100), A[1:25])
#'
#' permTestTx_results <- permTestTx(A, B, txdb, ntimes = 5)

permTestTx <- function(RS1 = NULL, RS2 = NULL, txdb = NULL, type = "mature", ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...) {
    RSL <- randomizeFeaturesTx(RS1, txdb, type, N = ntimes)
    permTestTx.results <- permTestTx_customAll(RSL = RSL, RS1 = RS1, RS2 = RS2, txdb = txdb, ev_function_1 = ev_function_1, ev_function_2= ev_function_2, pval_z = pval_z,...)
    return(permTestTx.results)
}

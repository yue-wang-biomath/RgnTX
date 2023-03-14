#' Perform permutation test
#' @export permTestTx_customAll
#' @importFrom stats sd
#' @importFrom stats pnorm
#'
#' @description Perform permutation test for evaluating spatial association between region sets. This permutation test function receives two region sets and a set of randomized region sets of one of them. It evaluates if there is an association between these two region sets.
#'
#' @usage permTestTx_customAll(RSL = NULL, RS1 = NULL, RS2 = NULL,
#' ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_t = FALSE, ...)
#'
#' @param RSL Randomized region sets of \code{RS1}. It should be a list object and each element should be in the \code{GRanges} or \code{GRangesList} format.
#' @param RS1 The region set. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param RS2 The region set to be compared with. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between RS1 and RS2. Default is  \code{overlapCountsTx}.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in RSL and RS2. Default is  \code{overlapCountsTx}.
#' @param pval_t Boolean. Default is FALSE. If FALSE, p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, p-value is calculated based on a t-test.
#' @param ... Any additional parameters needed.
#'
#' @details \code{permTestTx_customAll} will use evaluation function \code{ev_function_1} to calculate the test statistic between \code{RS1} and \code{RS2}, and use \code{ev_function_2} to evaluate the statistic between \code{RSL} and \code{RS2}. It will also return a p-value and a z-score.
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
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids1<- c("170")
#' RS1 <- randomizeTx(txdb = txdb, trans_ids = trans.ids1,
#'                     random_num = 20, random_length = 100)
#' RS2 <- randomizeTx(txdb = txdb, trans_ids = trans.ids1,
#'                     random_num = 20, random_length = 100)
#' trans.ids2 <-  c("170", "782", "974", "1364", "1387")
#' RSL <- randomizeTx(txdb = txdb, trans_ids = trans.ids2,
#'                     random_num = 20, random_length = 100, N = 10)
#' permTestTx_results <- permTestTx_customAll(RSL = RSL, RS1 = RS1, RS2 = RS2)
permTestTx_customAll <- function(RSL = NULL, RS1 = NULL, RS2 = NULL, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_t = FALSE, ...) {
    orig.ev <- ev_function_1(RS1, RS2,...)
    rand.ev <- lapply(RSL, function(x) {return(ev_function_2(RS2, x,...))})

    pval_zscore <- getPvalZscore(orig.ev, unlist(rand.ev), pval_t = pval_t)
    permTestTx.results <- list(RSL, RS1, RS2, orig.ev, unlist(rand.ev), pval_zscore[1], pval_zscore[2], length(unlist(rand.ev)))
    names(permTestTx.results) <- c("RSL", "RS1", "RS2", "orig.ev", "rand.ev", "pval", "zscore", "ntimes")
    class(permTestTx.results) <- "permTestTx.results"
    return(permTestTx.results)
}

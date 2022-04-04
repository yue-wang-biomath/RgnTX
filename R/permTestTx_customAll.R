#' Perform permutation test.
#' @export permTestTx_customAll
#' @importFrom stats sd
#' @importFrom stats pnorm
#'
#' @description Perform permutation test for evaluating spatial association between region sets. This permutation test function receives two region sets and a set of randomized region sets of one of them. It evaluates if there is an association between these two region sets.
#'
#' @usage permTestTx_customAll(RSL = NULL, RS1 = NULL, RS2 = NULL,
#' ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...)
#'
#' @param RSL Randomized region sets of \code{RS1}. It should be a list object with each element being in the \code{GRanges} or \code{GRangesList} format.
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param RS2 The region set to be compared with. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between RS1 and RS2. Default is overlapCountsTx.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in RSL and RS2. Default is overlapCountsTx.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... any additional parameters needed.
#'
#' @details \code{permTestTx_customAll} will use evaluation function to calculate the test statistic value between \code{RS1} and \code{RS2}, and the statistic values between \code{RSL} and \code{RS2}. It will also return a p-value and a z-score.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following items:
#' \itemize{
#' \item \bold{\code{RSL}}
#' \item \bold{\code{RS1}}
#' \item \bold{\code{RS2}}
#' \item \bold{\code{orig.ev}} The value of the test statistc between RS1 and RS2.
#' \item \bold{\code{rand.ev}} The values of the test statistc between each element in RSL and RS2.
#' \item \bold{\code{pval}} The p-value of the test.
#' \item \bold{\code{zscore}} The standard score of the test.
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
permTestTx_customAll <- function(RSL = NULL, RS1 = NULL, RS2 = NULL, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...) {
    orig.ev <- ev_function_1(RS1, RS2,...)
    rand.ev <- lapply(RSL, function(x) {
                return(ev_function_2(RS2, x,...))
    })
    rand.ev <- unlist(rand.ev)

    ntimes <- length(rand.ev)
    # orig.ev < rand.ev
    if (orig.ev < mean(rand.ev)) {
        zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE)) / sd(rand.ev, na.rm = TRUE), 4)
        xcoords <- rand.ev
        rand.mean <- mean(xcoords, na.rm = TRUE)
        rand.sd <- sd(xcoords, na.rm = TRUE)
        ntimes <- length(RSL)

        if (pval_z == FALSE){
            pval <- (sum(orig.ev >= rand.ev, na.rm = TRUE) + 1) / (ntimes + 1)
        }else{
            pval <- pnorm(orig.ev, mean(rand.ev), sd(rand.ev),lower.tail= FALSE)
        }
    }
    # orig.ev >= rand.ev
    if (orig.ev >= mean(rand.ev)) {
        zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE)) / sd(rand.ev, na.rm = TRUE), 4)
        xcoords <- rand.ev
        rand.mean <- mean(xcoords, na.rm = TRUE)
        rand.sd <- sd(xcoords, na.rm = TRUE)
        ntimes <- length(RSL)
        if (pval_z == FALSE){
            pval <- (sum(orig.ev <= rand.ev, na.rm = TRUE) + 1) / (ntimes + 1)
        }else{
            pval <- pnorm(orig.ev, mean(rand.ev), sd(rand.ev),lower.tail= TRUE)
        }
    }
    pval <- round(pval, 6)

    permTestTx.results <- list(
        RSL, RS1, RS2, orig.ev, rand.ev, pval, zscore,
        ntimes
    )
    names(permTestTx.results) <- c(
        "RSL", "RS1", "RS2", "orig.ev", "rand.ev",
        "pval", "zscore", "ntimes"
    )
    class(permTestTx.results) <- "permTestTx.results"
    return(permTestTx.results)
}

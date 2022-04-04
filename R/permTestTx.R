#' Perform permutation test.
#' @export permTestTx
#'
#' @description Perform permutation test for evaluating spatial association between features and a region set.
#'
#' @usage permTestTx(RS1 = NULL, RS2 = NULL, txdb = NULL, type = 'mature',
#' ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx,
#' pval_z = FALSE, ...)
#'
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param RS2 The region set to be compared with. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A txdb object.
#' @param type This argument receives options 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR'. It decides which types of space that the features being randomized into.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between RS1 and RS2. Default is overlapCountsTx.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in RSL and RS2. Default is overlapCountsTx.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... any additional parameters needed.
#'
#' @details \code{permTestTxIA} only needs users to input two region sets. It will automatically randomize the first region set into transcriptome.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following items:
#' \itemize{
#' \item \bold{\code{RSL}} Randomized region sets of \code{RS1}.
#' \item \bold{\code{RS1}} The region set to be randomized.
#' \item \bold{\code{RS2}} The region set to be compared with.
#' \item \bold{\code{orig.ev}} The value of the test statistic between RS1 and RS2.
#' \item \bold{\code{rand.ev}} The values of the test statistic between each element in RSL and RS2.
#' \item \bold{\code{pval}} The p-value of the test.
#' \item \bold{\code{zscore}} The standard score of the test.
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
    orig.ev <- ev_function_1(RS1, RS2, ...)
    RSL <- randomizeFeaturesTx(RS1, txdb, type, N = ntimes)

    rand.ev <- lapply(RSL, function(x) {
        return(ev_function_2(RS2, x, ...))
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

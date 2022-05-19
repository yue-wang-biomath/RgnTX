#' Perform permutation test
#' @export permTestTxIA_customPick
#'
#' @description Perform permutation test for evaluating spatial association between RNA features and a specified kind of regions. The latter is defined by the \code{customPick_function} argument input by users. The difference between this function and  \code{\link{permTestTx_customPick}} is that it is for RNA-related genomic features that have isoform ambiguity, i.e., features that one does not know which transcript they comes from.
#'
#' @usage permTestTxIA_customPick(RS1 = NULL, txdb = NULL, type = 'mature',
#' customPick_function = NULL, ntimes = 50,
#' ev_function_1 = overlapCountsTxIA, ev_function_2 = overlapCountsTx,  pval_z = FALSE, ...)
#'
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @param customPick_function A custom function needs to be inputted by users. The customPick function should have two arguments: a TxDb object and a character object of transcript ids. It returns a specified region over each transcript.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between \code{RS1} and \code{RS2}. Default is \code{overlapCountsTxIA}.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in \code{RSL} and \code{RS2}. Default is \code{overlapCountsTx}.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... Any additional parameters needed.
#'
#' @details \code{permTestTxIA_customPick} will assess the test statistic between \code{RS1} and each region in \code{RSL}, and the relation between \code{RS1} and \code{RS2}.
#' Each RNA feature is only mapped with a part of region on its transcript (picked by the \code{customPick_function}). The output \code{orig.ev} is the weighted counts between \code{RS1} and \code{RS2}. Each feature in \code{RS1} related to \code{n1} isoforms in \code{RS2} and overlapped with \code{n2} \code{RS2} regions will contribute a value of \code{n2/n1} to the total number of overlaps.
#' This test function also randomizes input features per transcript. The set of randomized results is outputted as \code{RSL}. The overlapping counts between each set in  \code{RSL} with \code{RS2} is outputted as \code{rand.ev}.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following information:
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
#'
#' permTestTx_results <- permTestTxIA_customPick(RS1 = RS1,
#'                                             txdb = txdb,
#'                                             type = 'mature',
#'                                             customPick_function = getStopCodon,
#'                                             ntimes = 5)
permTestTxIA_customPick <- function(RS1 = NULL, txdb = NULL, type = "mature",
                                    customPick_function = NULL, ntimes = 50, ev_function_1 = overlapCountsTxIA, ev_function_2 =overlapCountsTx, pval_z = FALSE, ...) {
    # RS1 should be a DNA feature in the GRanges format.
    trans.info <- getTransInfo(RS1, txdb)
    trans.id <- trans.info[, "trans_ID"]

    # Use input custom_function function to generate RS2.
    RS2 <- customPick_function(trans_ids = trans.id, txdb = txdb, ...)

    # Generate a list of randomized regions of RS1.
    RSL <- randomizeFeaturesTxIA(RS1, txdb, type, N = ntimes)

    # Evaluate the observed value.
    orig.ev <- ev_function_1(RS1, RS2,...)
    rand.ev <- lapply(RSL, function(x) {return(ev_function_2(x, RS2, ...))})
    pval_zscore <- getPvalZscore(orig.ev, unlist(rand.ev))
    permTestTx.results <- list(RSL, RS1, RS2, orig.ev, unlist(rand.ev), pval_zscore[1], pval_zscore[2], length(unlist(rand.ev)))
    names(permTestTx.results) <- c("RSL", "RS1", "RS2", "orig.ev", "rand.ev", "pval", "zscore", "ntimes")
    class(permTestTx.results) <- "permTestTx.results"
    return(permTestTx.results)

}

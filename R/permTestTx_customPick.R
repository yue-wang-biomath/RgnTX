#' Perform permutation test
#' @export permTestTx_customPick
#'
#' @description Perform permutation test for evaluating spatial association between a feature set and the customPick regions. The latter is defined by the \code{customPick_function} argument provided by users.
#'
#' @usage permTestTx_customPick(RS1 = NULL, txdb = NULL, type = "mature",
#' customPick_function = NULL, ntimes = 50, ev_function_1 = overlapCountsTx,
#' ev_function_2 = overlapCountsTx,  pval_z = FALSE, ...)
#'
#' @param RS1 The feature set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @param customPick_function A custom function needs to be inputted by users. The custom function should have two arguments: a TxDb object and a character object of transcript ids. It returns a part of region of each transcript.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between \code{RS1} and \code{RS2}. Default is \code{overlapCountsTx}.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in \code{RSL} and \code{RS2}. Default is \code{overlapCountsTx}.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... Any additional parameters needed.
#'
#' @details
#' Each feature in \code{RS1} is only mapped with the customPick regions over its transcript (picked by the \code{customPick_function}). The output \code{orig.ev} is the number of features that have overlap with its customPick region.
#' The set of randomized region sets is outputted as \code{RSL}. The overlapping counts between each set in \code{RSL} with \code{RS2} is outputted as \code{rand.ev}.
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
#' trans.ids <- sample(names(exons.tx0), 100)
#' RS1 <- randomizeTx(txdb, trans.ids, random_num = 100,
#' random_length = 200, type = 'CDS')
#' getCDS = function(txdb, trans.id){
#' cds.tx0 <- cdsBy(txdb, use.names=FALSE)
#'     cds.names <- as.character(intersect(names(cds.tx0), trans.id))
#'     cds = cds.tx0[cds.names]
#'     return(cds)
#' }
#'
#' permTestTx_results <- permTestTx_customPick(RS1,txdb,
#' customPick_function = getCDS, ntimes = 5)

permTestTx_customPick <- function(RS1 = NULL, txdb = NULL, type = "mature", customPick_function = NULL, ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...) {
    if (is(RS1, "GRanges")) {
        trans.id <- RS1$transcriptsHits
    }
    if (is(RS1, "CompressedGRangesList")) {
        trans.id <- names(RS1)
    }
    # Use input custom_function function to generate RS2.
    RS2 <- customPick_function(unique(trans.id),txdb = txdb, ...)
    orig.ev <- ev_function_1(RS1, RS2, ...)
    RSL <- randomizeFeaturesTx(RS1, txdb, type, N = ntimes)

    rand.ev <- lapply(RSL, function(x) {return(ev_function_2(x, RS2, ...))})
    pval_zscore <- getPvalZscore(orig.ev, unlist(rand.ev))
    permTestTx.results <- list(RSL, RS1, RS2, orig.ev, unlist(rand.ev), pval_zscore[1], pval_zscore[2], length(unlist(rand.ev)))
    names(permTestTx.results) <- c("RSL", "RS1", "RS2", "orig.ev", "rand.ev", "pval", "zscore", "ntimes")
    class(permTestTx.results) <- "permTestTx.results"
    return(permTestTx.results)
}

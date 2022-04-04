#' Perform permutation test.
#' @export permTestTxIA_customPick
#'
#' @description Perform permutation test for evaluating spatial association between RNA features and a specified kind of regions. The latter is defined by the \code{customPick_function} argument input by users. The difference between this function and  \code{\link{permTestTx_customPick}} is that it is for RNA-related genomic features that have isoform ambiguity, i.e., features that one does not know which transcript they comes from.
#'
#' @usage permTestTxIA_customPick(RS1 = NULL, txdb = NULL, type = 'mature',
#' customPick_function = NULL, ntimes = 50,
#' ev_function_1 = overlapCountsTxIA, ev_function_2 = overlapCountsTx,  pval_z = FALSE, ...)
#'
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A txdb object.
#' @param type This argument receives options 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR'. It decides which types of space that the features being randomized into.
#' @param customPick_function A custom function needs to be input by users. The custom function should have two arguments: a txdb object and a character object of transcript ids. It returns a specified region over each transcript.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between RS1 and RS2. Default is overlapCountsTxIA.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in RSL and RS2. Default is overlapCountsTx.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... any additional parameters needed.
#'
#' @details \code{permTestTxIA_customPick} will assess the relation between RS1 and RSL, and the relation between \code{RS1} and \code{RS2}.
#' Each RNA feature is only mapped with a part of region on its transcript (picked by the \code{customPick_function}). The output \code{orig.ev} is the weighted counts between \code{RS1} and \code{RS2}. For example, if a feature in \code{RS1} have overlap
#' This test function also randomizes input features per transcript. The set of randomized results is outputted as \code{RSL}. The overlapping counts between each set in  \code{RSL} with \code{RS2} is outputted as \code{rand.ev}.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following information:
#' \itemize{
#' \item \bold{\code{RSL}} Randomized region sets of \code{RS1}.
#' \item \bold{\code{RS1}} Features.
#' \item \bold{\code{RS2}} Region set to be compared with.
#' \item \bold{\code{orig.ev}} The value of weighted overlapping counts between RS1 and RS2.
#' \item \bold{\code{rand.ev}} The values of overlapping counts between each element in RSL and RS2.
#' \item \bold{\code{pval}} The p-value of the test.
#' \item \bold{\code{zscore}} The standard score of the test.
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

    # Evaluation step.
    orig.ev <- ev_function_1(RS1, RS2,...)

    # Generate a list of randomized regions of RS1.
    RSL <- randomizeFeaturesTxIA(RS1, txdb, type, N = ntimes)

    # Calculate overlapping counts between RSL and RS2.
    rand.ev <- lapply(RSL, function(x) {
        return(ev_function_2(x, RS2, ...))
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

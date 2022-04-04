#' Perform permutation test.
#' @export permTestTx_customPick
#'
#' @description Perform permutation test for evaluating spatial association between RNA features and a specified kind of regions. The latter is defined by the \code{customPick_function} argument input by users.
#'
#' @usage permTestTx_customPick(RS1 = NULL, txdb = NULL, type = "mature",
#' customPick_function = NULL, ntimes = 50, ev_function_1 = overlapCountsTx,
#' ev_function_2 = overlapCountsTx,  pval_z = FALSE, ...)
#'
#' @param RS1 The region set to be randomized. It should be in the \code{GRanges} or \code{GRangesList} format.
#' @param txdb A txdb object.
#' @param type This argument receives options 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR'. It decides which types of space that the features being randomized into.
#' @param customPick_function A custom function needs to be inputted by users. The custom function should have two arguments: a txdb object and a character object of transcript ids. It returns a part of region of each transcript.
#' @param ntimes Randomization times.
#' @param ev_function_1 Evaluation function defines what statistic to be tested between RS1 and RS2. Default is overlapCountsTx_one2one.
#' @param ev_function_2 Evaluation function defines what statistic to be tested between each element in RSL and RS2. Default is overlapCountsTx_one2one.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#' @param ... any additional parameters needed.
#'
#' @details
#' The transcript where each RNA feature comes from should be available. Each RNA feature is only mapped with a part of region on its transcript (picked by the \code{customPick_function}). The output \code{orig.ev} is the number of features that have overlap with its specified region.
#' This test function also randomizes input features per transcript. The set of randomized results is outputted as \code{RSL}. The overlapping counts between each set in \code{RSL} with \code{RS2} is outputted as \code{rand.ev}.
#'
#' @return
#' A list object, which is defined to be \code{permTestTx.results} class. It contains the following information:
#' \itemize{
#' \item \bold{\code{RSL}} Randomized region sets of \code{RS1}.
#' \item \bold{\code{RS1}} Features.
#' \item \bold{\code{RS2}} Region set to be compared with.
#' \item \bold{\code{orig.ev}} The value of overlapping counts between RS1 and RS2.
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
#' exons.tx0 <- exonsBy(txdb)
#' trans.ids <- sample(names(exons.tx0), 100)
#' RS1 <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 200, type = 'CDS')
#' getCDS = function(txdb, trans.id){
#' cds.tx0 <- cdsBy(txdb, use.names=FALSE)
#'     cds.names <- as.character(intersect(names(cds.tx0), trans.id))
#'     cds = cds.tx0[cds.names]
#'     return(cds)
#' }
#'
#' permTestTx_results <- permTestTx_customPick(RS1,txdb,
#'                                             customPick_function = getCDS,
#'                                             ntimes = 5)

permTestTx_customPick <- function(RS1 = NULL, txdb = NULL, type = "mature", customPick_function = NULL, ntimes = 50, ev_function_1 = overlapCountsTx, ev_function_2 = overlapCountsTx, pval_z = FALSE, ...) {
    if (is(RS1, "GRanges")) {
        trans.id <- RS1$transcriptsHits
    }
    if (is(RS1, "CompressedGRangesList")) {
        trans.id <- names(RS1)
    }

    # Use input custom_function function to generate RS2.
    RS2 <- customPick_function(unique(trans.id),txdb = txdb, ...)

    # Evaluation step.
    orig.ev <- ev_function_1(RS1, RS2, ...)

    # Generate a list of randomized regions of RS1.
    RSL <- randomizeFeaturesTx(RS1, txdb, type, N = ntimes)

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

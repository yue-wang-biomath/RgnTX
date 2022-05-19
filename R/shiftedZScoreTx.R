#' Calculate shifted z scores
#' @export shiftedZScoreTx
#'
#' @description Calculate shifted z scores for permutation test results.
#'
#' @usage shiftedZScoreTx(permTestTx_results = NULL, txdb = NULL,
#' window = 200, step = 20, ev_function_1 = overlapCountsTx, ...)
#'
#' @param permTestTx_results A \code{permTestTx.results} object.
#' @param txdb A TxDb object.
#' @param window The window of the whole shifting.
#' @param step The step of each shifting.
#' @param ev_function_1 Evaluation function. Default is \code{overlapCountsTx}.
#' @param ... Any additional parameters needed.
#'
#' @return
#' A list object, which is defined to be \code{shitedZScore.results} class. It contains the following items:
#' \itemize{
#' \item \bold{\code{shifted.z.scores:}} Standard z-scores after shifting.
#' \item \bold{\code{window:}} Window of the whole shifting.
#' \item \bold{\code{step:}} Step of each shifting.
#' \item \bold{\code{original.z.score:}} Original standard score.
#' }
#'
#' @seealso \code{\link{plotShiftedZScoreTx}}
#'
#' @details see examples in \code{\link{plotShiftedZScoreTx}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' file <- system.file(package="RgnTX", "extdata/m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' RS1 <- m6A_sites_data[1:500]
#' permTestTx_results <- permTestTxIA_customPick(RS1 = RS1,
#'                                         txdb = txdb,
#'                                         customPick_function = getStopCodon,
#'                                         ntimes = 5)
#' shiftedZScoreTx_results <- shiftedZScoreTx(permTestTx_results, txdb = txdb,
#'                                         window = 2000,
#'                                         step = 200,
#'                                         ev_function_1 = overlapCountsTxIA)
shiftedZScoreTx <- function(permTestTx_results = NULL, txdb = NULL, window = 200, step = 20, ev_function_1 = overlapCountsTx, ...) {
    if(!is(permTestTx_results, 'permTestTx.results')){
        stop("Argument permTestTx_results must be a permTestTx.results object.")
    }
    RS1 <- permTestTx_results$RS1
    RS2 <- permTestTx_results$RS2
    rand.ev <- permTestTx_results$rand.ev
    original.z.score <- permTestTx_results$zscore

    A <- RS1
    B <- RS2
    mean.permuted <- mean(rand.ev)
    sd.permuted <- sd(rand.ev)

    if(length(window) == 0 || length(step) == 0){
        window <- 5 * mean(sum(width(A)))
        step <- floor(window / 10)
    }
    num.steps <- floor(window / step)

    shifts <- seq_len(num.steps) * step
    shifts <- c(rev(-1 * shifts), 0, shifts)

    if (is(A, "CompressedGRangesList")) {
        A.unlist <- GRangesList2GRanges(A)
    } else {
        A.unlist <- A
    }

    shifted.z.score <- function(shift) {
        shifted.A <- shift(A.unlist, shift)
        shifted.evaluation <- ev_function_1(shifted.A, B, ...)
        shifted.z.score <- (shifted.evaluation - mean.permuted) / sd.permuted
        return(shifted.z.score)
    }

    shifted <- lapply(shifts, shifted.z.score)
    shifted.z.scores <- do.call(c, shifted)
    shifted.z.scores[num.steps+1] <- original.z.score

    shitedZScoresTx_results <- list(shifted.z.scores, shifts, window, original.z.score)
    names(shitedZScoresTx_results) <- c(
        "shifted.z.scores", "shifts", "window",
        "original.z.score"
    )

    class(shitedZScoresTx_results) <- "shitedZScoreTx.results"
    return(shitedZScoresTx_results)
}

#' Calculate shifted z scores over transcripts
#' @export shiftedZScoreTx
#'
#' @description Calculate shifted z scores for permutation test results. It provides mode of shifting regions of interest over mRNA space.
#'
#' @usage shiftedZScoreTx(permTestTx_results = NULL, txdb = NULL, type = "normal",
#' window = 200, step = 20, ev_function_1 = overlapCountsTx, ...)
#'
#' @param permTestTx_results A \code{permTestTx.results} object.
#' @param txdb A TxDb object.
#' @param type If type is "mature", regions of interest will be shifted over the mature mRNA space according to the input txdb object. Otherwise ROS will be shifted over genomic space.
#' @param window The window of the whole shifting.
#' @param step The step of each shifting.
#' @param ev_function_1 Evaluation function. Default is \code{overlapCountsTx}.
#' @param ... Any additional parameters needed.
#'
#' @return
#' A list object, which is defined to be \code{shiftedZScoreTx.results} class. It contains the following items:
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
#' txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' file <- system.file(package="RgnTX", "extdata", "m6A_sites_data.rds")
#' m6A_sites_data <- readRDS(file)
#' RS1 <- m6A_sites_data[1:100]
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' permTestTx_results <- permTestTxIA_customPick(RS1 = RS1,
#'                                               txdb = txdb,
#'                                               type = "mature",
#'                                               customPick_function = getStopCodon,
#'                                               ntimes = 1)
#' shiftedZScoreTx_results <- shiftedZScoreTx(permTestTx_results,txdb,
#'                                            type = 'mature',
#'                                            window = 500,
#'                                            step = 50,
#'                                            ev_function_1 = overlapCountsTxIA)
shiftedZScoreTx <- function(permTestTx_results = NULL,
                            txdb = NULL,
                            type = "normal",
                            window = 200,
                            step = 20,
                            ev_function_1 = overlapCountsTx, ...) {
    if(!is(permTestTx_results, "permTestTx.results")){
        stop("Argument permTestTx_results must be a permTestTx.results object.")
    }
    RS1 <- permTestTx_results$RS1
    RS2 <- permTestTx_results$RS2
    rand.ev <- permTestTx_results$rand.ev
    original.z.score <- permTestTx_results$zscore

    A <- RS1
    B <- RS2
    B = adjustGroup(B)
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

    if(type == 'mature'){
        dataframe.B<- data.frame(B)
        strand.B<- dataframe.B[, 'strand']
        strand.B <- strand.B[unlist(lapply(unique(dataframe.B[,1]), function(x){return(min(which(dataframe.B[, 1] == x)))}))]
        width.B <- sum(width(B))
        start.B <- 1:length(B)
        start.B[strand.B == '+'] <- min(start(B[strand.B == '+']))
        start.B[strand.B == '-'] <- max(end(B[strand.B == '-']))
        exons.tx0 <- exonsBy(txdb)
        exons.tx.B <- exons.tx0[names(B)]

        shiftTxB1 <- function(distance){
            distance <- replicate(length(B), distance)
            # calculate start.shift.B
            shift.regions.temp <- shiftExonTx(exons.tx.B, start.B, distance)
            start.shift.B <- 1:length(B)
            if(distance[1] > 0){
                start.shift.B[strand.B == '+'] <- max(end(shift.regions.temp[strand.B == '+']))
                start.shift.B[strand.B == '-'] <- min(start(shift.regions.temp[strand.B == '-']))
                shift.regions <- shiftExonTx(exons.tx.B, start.shift.B, width.B)
            }
            if(distance[1] < 0){
                start.shift.B[strand.B == '+'] <- min(start(shift.regions.temp[strand.B == '+']))
                start.shift.B[strand.B == '-'] <- max(end(shift.regions.temp[strand.B == '-']))
                shift.regions <- shiftExonTx(exons.tx.B, start.shift.B, -width.B)
            }
            return(shift.regions)
        }
        shiftTxB <- shiftTxB1
    } else {
        shiftTxB2 <- function(width) {
            shift.regions <- shift(unlist(B), width)
            return(shift.regions)
        }
        shiftTxB <- shiftTxB2
    }

    # Shift regions B
    shifts_left <- shifts[1:num.steps]
    shifts_right <- shifts[(num.steps+2):length(shifts)]
    shifted.B_left <- lapply(shifts_left, shiftTxB)
    shifted.B_right <- lapply(shifts_right, shiftTxB)
    shifted.B <- c(shifted.B_left,B, shifted.B_right)
    # Evaluate association between shifted B with A
    ev_function_1_fun <- function(B){
        return(ev_function_1(A, B))
    }
    # Evaluation results
    shifted.evaluation <- unlist(lapply(shifted.B, ev_function_1_fun))
    shifted.z.scores <- (shifted.evaluation - mean.permuted) / sd.permuted

    shiftedZScoreTx_results <- list(shifted.z.scores, shifts, window, original.z.score)
    names(shiftedZScoreTx_results) <- c(
        "shifted.z.scores", "shifts", "window",
        "original.z.score"
    )

    class(shiftedZScoreTx_results) <- "shiftedZScoreTx.results"
    return(shiftedZScoreTx_results)
}

#' Plot shifted z scores
#' @export plotShiftedZScoreTx
#'
#' @description Plot shifted z scores for permutation test results.
#'
#' @usage plotShiftedZScoreTx(shitedZScoresTx_results)
#'
#' @param shitedZScoresTx_results A \code{shitedZScoreTx.results} object.
#'
#' @return A plot.
#'
#' @seealso \code{\link{shiftedZScoreTx}}
#' @examples
#' file <- system.file(package="RgnTX", "extdata", "shiftedZScoreTx_results1.rds")
#' shiftedZScoreTx_results <- readRDS(file)
#' p1 <- plotShiftedZScoreTx(shiftedZScoreTx_results)
#' p1


plotShiftedZScoreTx <- function(shitedZScoresTx_results) {
    if(!is(shitedZScoresTx_results, 'shitedZScoreTx.results')){
        stop("Argument shitedZScoresTx_results must be a shitedZScoreTx.results object.")
    }
    shifted.z.scores <- shitedZScoresTx_results$shifted.z.scores
    shifts <- shitedZScoresTx_results$shifts
    window <- shitedZScoresTx_results$window
    original.z.score <- shitedZScoresTx_results$original.z.score
    data <- data.frame(shifts = shifts, shifted.z.scores = shifted.z.scores)

    ymax <- max(shifted.z.scores, 2)
    ymin <- min(shifted.z.scores, -2)

    p1 <- ggplot(data, aes(x = shifts, y = shifted.z.scores)) +
        theme(
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 1),
            panel.grid.major = element_line(colour = "grey", linetype = 9, size = 0.2),
            axis.ticks = element_blank(),
            line = element_line(colour = "white", size = 0.5, linetype = 9, lineend = "butt"),
            legend.position = "bottom",
            axis.text = element_text(size = 11),
            legend.text = element_text(size = 11),
            axis.title.y = element_text(size = 11, hjust = 0.6),
            title = element_text(size = 11, face = "bold")
        ) +
        geom_line() +
        ylim(ymin, ymax)
    return(p1)
}

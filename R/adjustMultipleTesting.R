#' Calculate adjusted pvalues
#' @export adjustMultipleTesting
#' @importFrom stats p.adjust
#'
#' @description Calculate adjusted pvalues of multiple testing.
#'
#' @usage adjustMultipleTesting(pval_list, alpha)
#'
#' @param pval_list Pvalues from multiple testing
#' @param alpha threshold
#'
#' @return
#' A list object, which contains the following two items:
#' \itemize{
#' \item \bold{\code{Table}} A data frame object with columns of rank, pvals, adjusted_pvals and whether to reject null hypothesis.
#' \item \bold{\code{Proportion}} The proportion of the tests that reject the null hypothesis H0 among all the tests.
#' }
#'
#' @examples
#' file <- system.file(package="RgnTX", "extdata/multi_pvals.rds")
#' multi_pvals <- readRDS(file)
#' adjustMultipleTesting(multi_pvals[, 1], 0.05)

adjustMultipleTesting = function(pval_list, alpha= 0.05){
    pval_N <- length(pval_list)
    pval_list <- sort(pval_list)
    pval_list_adjust <- p.adjust(pval_list, "fdr")
    index_reject <- max(which(pval_list_adjust <= alpha))
    proportion_reject <- index_reject/pval_N
    reject_H0 <- c(replicate(index_reject, 'Yes'), replicate(pval_N - index_reject, 'No'))
    data_frame_p <- data.frame(rank = 1:pval_N,
                               pval = pval_list,
                               adjusted_pval = pval_list_adjust,
                               reject_H0 = reject_H0)
    adjusted_results <- list(data_frame_p, proportion_reject)
    names(adjusted_results) <- c('Table','Proportion')
    return(adjusted_results)
    }

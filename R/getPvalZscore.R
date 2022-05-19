#' getPvalZscore
#' @export getPvalZscore
#'
#' @description Calculate a p-value and z-score based on observed value and random evaluation values.
#'
#' @usage getPvalZscore(orig.ev, rand.ev, pval_z = FALSE)
#'
#' @param orig.ev Observed value.
#' @param rand.ev Random evaluation values.
#' @param pval_z Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a z-test.
#'
#' @return A p-value and a z-score.
#'
#' @examples
#' orig.ev <- 5
#' random.ev <- rnorm(10, mean =2, sd = 1)
#' getPvalZscore(orig.ev, random.ev)
getPvalZscore <- function(orig.ev, rand.ev, pval_z = FALSE){
    # orig.ev < rand.ev
    if (orig.ev < mean(rand.ev)) {
        zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE)) / sd(rand.ev, na.rm = TRUE), 4)
        xcoords <- rand.ev
        rand.mean <- mean(xcoords, na.rm = TRUE)
        rand.sd <- sd(xcoords, na.rm = TRUE)
        ntimes <- length(rand.ev)

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
        ntimes <- length(rand.ev)
        if (pval_z == FALSE){
            pval <- (sum(orig.ev <= rand.ev, na.rm = TRUE) + 1) / (ntimes + 1)
        }else{
            pval <- pnorm(orig.ev, mean(rand.ev), sd(rand.ev),lower.tail= TRUE)
        }
    }
    pval <- round(pval, 6)
    return(c(pval, zscore))
}

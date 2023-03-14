################################################################################
# UTILITY FUNCTIONS
################################################################################
# These functions remain interior to the package (not exported)

#' getFormatCorrect
#' @description This function makes sure the two input region sets are in the correct format required by RgnTX evaluation functions.
#' @usage getFormatCorrect(A, B)
#' @param A Region set 1. A Granges or GRangesList object.
#' @param B Region set 2. A Granges or GRangesList object.
#' @return A \code{list} object.
getFormatCorrect <- function(A, B){
    if (is(A, "GRangesList")) {
        A <- GRangesList2GRanges(A)
    }
    if (is(B, "GRangesList")) {
        B <- GRangesList2GRanges(B)
    }

    if (!is(A, "GRanges")) {
        stop("A should be either GRanges or GRangesList.")
    }
    if (!is(B, "GRanges")) {
        stop("B should be either GRanges or GRangesList.")
    }

    if (length(A$group) == 0) {
        A$group <- seq_len(length(A))
    }

    if (length(B$group) == 0) {
        B$group <- seq_len(length(B))
    }
    return(list(A,B))
}


#' Get permutation space by specifying transcript ids
#' @description This function returns 5'UTR/CDS/3'UTR/mRNA/full part of transcriptome regions grouped by corresponding transcript ids.
#' @usage getPermSpaceByTxID(trans_ids = "all", txdb, type = "mature")
#' @param trans_ids A character object. The transcript ids. Default is "all". If it takes the default value "all", the space that users get will be the whole transcriptome.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @return A \code{GRangesList} object.
#' @seealso \code{\link{getPermSpaceByType}}, \code{\link{getPermSpaceByFeatures}}
getPermSpaceByTxID <- function(trans_ids = "all", txdb, type = "mature") {

    # make sure the inputs are reasonable
    trans.ids <- trans_ids
    if (!is.character(trans.ids)) {
        stop("trans.ids must be character.")
    }

    regions.A <- getPermSpaceByType(txdb, type = type)

    if (trans.ids[1] == "all") {
        trans.ids <- names(regions.A)
    } else {
        trans.ids <- as.character(trans.ids)
        trans.ids <- trans.ids[is.element(trans.ids, names(regions.A))]
    }

    if (length(trans.ids) != 0) {
        regions.A <- regions.A[trans.ids]
    } else {
        regions.A <- NULL
    }

    return(regions.A)
}


#' Get permutation space for features
#' @description This function returns a default permutation space for features with isoform ambiguity. The default permutation space of a feature is the aggregate of the multiple transcripts it may overlap with. It requires the input feature to be \code{GRanges} format.
#' @usage getPermSpaceByFeatures(features, txdb, type = "mature")
#' @param features A \code{GRanges} object.
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of regions over transcriptome.
#' @return
#' A list object, which contains two elements.
#' \itemize{
#' \item \bold{\code{perm.space:}} A \code{GRangesList} object that includes all the transcripts input features may overlap with.
#' \item \bold{\code{index:}} It contains a series of numbers indicating which feature these transcripts are respectively associated with.
#' }
#' @seealso \code{\link{getPermSpaceByTxID}}, \code{\link{getPermSpaceByType}}
getPermSpaceByFeatures <- function(features, txdb, type = "mature") {

    # make sure the inputs are reasonable

    if (!is.character(type)) {
        stop("type must be character.")
    }

    regions.A <- getPermSpaceByType(txdb, type = type)
    A <- features

    if (!is(A, "GRanges")) {
        stop("features must be a GRanges object.")
    } else {
        trans_info <- getTransInfo(A, txdb)

        trans_id <- as.character(trans_info[, "trans_ID"])
        index_methyl <- trans_info[, "index_features"]
        index_trans <- trans_info[, "index_trans"]

        regions.A <- regions.A[trans_id]
    }
    return.list <- list(regions.A, index_methyl)
    names(return.list) <- c("perm.space", "index")
    return(return.list)
}


#' Get permutation space by specifying type
#' @importFrom GenomicFeatures fiveUTRsByTranscript threeUTRsByTranscript cdsBy transcriptsBy exonsBy
#' @description This function can return 5'UTR/CDS/3'UTR/mRNA/full part of transcriptome regions, following the format required by the main permutation test functions.
#' @usage getPermSpaceByType(txdb, type = "mature")
#' @param txdb A TxDb object.
#' @param type A character object. Default is "mature". It accepts options "mature", "full", "fiveUTR", "CDS" or "threeUTR", with which one can get corresponding types of transcriptome regions.
#' @return A \code{GRangesList} object.
#' @seealso \code{\link{getPermSpaceByTxID}}, \code{\link{getPermSpaceByFeatures}}
getPermSpaceByType <- function(txdb, type = "mature") {

    # make sure the inputs are reasonable
    if (!is.character(type)) {
        stop("type must be character.")
    }

    regions.A <- NULL
    if (type == "mature") {
        exons.tx0 <- exonsBy(txdb, by = c("tx", "gene"), use.names = FALSE)
        regions.A <- exons.tx0
    }
    if (type == "fiveUTR") {
        fiveUTR.tx0 <- fiveUTRsByTranscript(txdb, use.names = FALSE)
        regions.A <- fiveUTR.tx0
    }
    if (type == "threeUTR") {
        threeUTR.tx0 <- threeUTRsByTranscript(txdb, use.names = FALSE)
        regions.A <- threeUTR.tx0
    }
    if (type == "CDS") {
        cds.tx0 <- cdsBy(txdb, use.names = FALSE)
        regions.A <- cds.tx0
    }
    if (type == "full") {
        trans.tx0 <- transcriptsBy(txdb, "gene")
        trans.tx0.df <- data.frame(trans.tx0)
        trans.tx0.id <- trans.tx0.df[, "tx_id"]
        # RefSeqID
        RefSeqID <- as.character(trans.tx0.id)
        # targetName
        targetName <- trans.tx0.df[, "seqnames"]
        targetName <- as.character(targetName)
        # strand
        strand <- trans.tx0.df[, "strand"]
        strand <- as.character(strand)
        # blockSizes
        blockSizes <- trans.tx0.df[, "width"]

        # targetStart
        targetStart <- trans.tx0.df[, "start"]

        regions.A <- vector2GRangesList(RefSeqID, targetName, strand, blockSizes, targetStart)
    }

    if (is.null(regions.A)) {
        stop("type must come from one of these options: 'mature', 'full', 'fiveUTR', 'CDS' or 'threeUTR'.")
    }

    return(regions.A)
}


#' getPvalZscore
#' @description Calculate a p-value and z-score based on observed value and random evaluation values.
#' @usage getPvalZscore(orig.ev, rand.ev, pval_t = FALSE)
#' @param orig.ev Observed value.
#' @param rand.ev Random evaluation values.
#' @param pval_t Boolean. Default is FALSE. If FALSE, the p-value is calculated based on the number of random evaluations is larger or less than the initial evaluation. If TRUE, the p-value is calculated based on a t-test.
#' @return A p-value and a z-score.
#' @importFrom stats t.test
getPvalZscore <- function(orig.ev, rand.ev, pval_t = FALSE){
    # orig.ev < rand.ev
    if (orig.ev < mean(rand.ev)) {
        zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE)) / sd(rand.ev, na.rm = TRUE), 4)
        xcoords <- rand.ev
        rand.mean <- mean(xcoords, na.rm = TRUE)
        rand.sd <- sd(xcoords, na.rm = TRUE)
        ntimes <- length(rand.ev)

        if (pval_t == FALSE){
            pval <- (sum(orig.ev >= rand.ev, na.rm = TRUE) + 1) / (ntimes + 1)
        }else{
            pval <- t.test(rand.ev - orig.ev, alternative = "greater")$p.value
        }
    }
    # orig.ev >= rand.ev
    if (orig.ev >= mean(rand.ev)) {
        zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE)) / sd(rand.ev, na.rm = TRUE), 4)
        xcoords <- rand.ev
        rand.mean <- mean(xcoords, na.rm = TRUE)
        rand.sd <- sd(xcoords, na.rm = TRUE)
        ntimes <- length(rand.ev)
        if (pval_t == FALSE){
            pval <- (sum(orig.ev <= rand.ev, na.rm = TRUE) + 1) / (ntimes + 1)
        }else{
            pval <- t.test(orig.ev - rand.ev, alternative = "greater")$p.value
        }
    }
    pval <- round(pval, 6)
    return(c(pval, zscore))
}

#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
adjustGroup = function(grl){
    dataframe.grl <- data.frame(grl)
    strand.grl <- dataframe.grl[, 'strand']
    strand.grl <- strand.grl[unlist(lapply(unique(dataframe.grl[,1]), function(x){return(min(which(dataframe.grl[, 1] == x)))}))]

    grl_start <- IntegerList(lapply(start(grl), function(x){return(sort(x))}))
    grl_start[strand.grl == '-'] <- IntegerList(lapply(grl_start[strand.grl == '-'], function(x){return(rev(sort(x)))}))
    grl_end <- IntegerList(lapply(end(grl), function(x){return(sort(x))}))
    grl_end[strand.grl == '-'] <- IntegerList(lapply(grl_end[strand.grl == '-'], function(x){return(rev(sort(x)))}))
    grl_start_temp <- IntegerList(lapply(start(grl), function(x){return(replicate(length(x), 0))}))
    start(grl) <- grl_start_temp
    end(grl) <- grl_end
    start(grl) <- grl_start
    return(grl)
}

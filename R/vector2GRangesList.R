#' vector2GRangesList
#' @export vector2GRangesList
#'
#' @description Generate \code{GRangesList} object from vectors. The output region set follows the format required by the main permutation test functions.
#'
#' @usage vector2GRangesList(RefSeqID, targetName, strand, blockSizes, targetStart)
#'
#' @param RefSeqID The name of each element.
#' @param targetName The seqnames of each range.
#' @param strand The strand of each range.
#' @param blockSizes The width of each range.
#' @param targetStart The stard coordinate of each range.
#'
#' @return A \code{GRangesList} object.
#'
#' @seealso \code{\link{GRanges2GRangesList}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.tx0 <- transcriptsBy(txdb, "gene")[1:2]
#' trans.tx0.df <- data.frame(trans.tx0)
#' trans.tx0.id <- trans.tx0.df[, "tx_id"]
#' # RefSeqID
#' RefSeqID <- as.character(trans.tx0.id)
#' # targetName
#' targetName <- trans.tx0.df[, "seqnames"]
#' targetName <- as.character(targetName)
#' # strand
#' strand <- trans.tx0.df[, "strand"]
#' strand <- as.character(strand)
#' # blockSizes
#' blockSizes <- trans.tx0.df[, "width"]
#' # targetStart
#' targetStart <- trans.tx0.df[, "start"]
#'
#' trans.list <- vector2GRangesList(RefSeqID, targetName, strand,
#' blockSizes, targetStart)

vector2GRangesList <- function(RefSeqID, targetName, strand, blockSizes, targetStart
){

targetStart <- as.character(targetStart)
targetStart <- lapply(targetStart, function(x) {
    process1 <- gsub("[(.*)]", "", x)
    process2 <- gsub("c", "", process1)
    process3 <- gsub(" ", "", process2)
    process4 <- gsub(":", ",", process3)
    process5 <- paste0(process4, ",")
    return(process5)
})
targetStart <- unlist(targetStart)

blockSizes <- as.character(blockSizes)
blockSizes <- lapply(blockSizes, function(x) {
    process1 <- gsub("[(.*)]", "", x)
    process2 <- gsub("c", "", process1)
    process3 <- gsub(" ", "", process2)
    process4 <- gsub(":", ",", process3)
    process5 <- paste0(process4, ",")
    return(process5)
})
blockSizes <- unlist(blockSizes)

A.frags <- cbind(RefSeqID, targetName, strand, targetStart, blockSizes)
A.frags <- data.frame(A.frags)
A.list <- with(A.frags, makeGRangesListFromFeatureFragments(
    seqnames = targetName,
    fragmentStarts = targetStart, fragmentWidths = blockSizes, strand = strand
))

if (length(which(is.na(RefSeqID) == FALSE)) != 0) {
    names(A.list) <- A.frags$RefSeqID
}

return(A.list)
}

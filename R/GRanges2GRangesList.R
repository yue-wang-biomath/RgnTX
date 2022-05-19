#' Convert a GRanges object to a GRangesList object
#' @export GRanges2GRangesList
#' @importFrom GenomicRanges GRanges
#'
#' @description Convert a \code{GRanges} object to a \code{GRangesList} object. The output region set follows the format required by the main permutation test functions.
#'
#' @usage GRanges2GRangesList(A = NULL)
#'
#' @param A A \code{GRanges} object.
#'
#' @details If input \code{GRanges} object has a metadata named as "group",
#' ranges having the same group number represent a region. If not, a range is a region.
#' A region in the input set will be outputted as a list element IN returned
#' \code{GRangesList} object.
#'
#' @return A \code{GRangesList} object.
#'
#' @seealso \code{\link{GRanges2GRangesList}}
#'
#' @examples
#' library(GenomicRanges)
#' GRanges.object <- GRanges(
#'     Rle(c("chr2", "chr2", "chr1", "chr3")),
#'     IRanges(1:4, width = 5)
#' )
#' # Assign the first and the second ranges to the same element.
#' GRanges.object$group <- c(1, 1, 2, 3)
#' GRangesList.object <- GRanges2GRangesList(GRanges.object)
GRanges2GRangesList <- function(A = NULL) {
    if (length(is.na(A$group)) == 0) {
        A$group <- seq_len(length(A))
    }

    metadata <- mcols(A)
    trans.id <- metadata$transcriptsHits
    group.name <- metadata$group
    group.start <- start(A)
    group.width <- width(A)
    group.index <- lapply(seq_len(max(group.name)), function(x) {
        return(which(group.name == x))
    })
    group.min <- unlist(lapply(seq_len(max(group.name)), function(x) {
        return(min(which(group.name == x)))
    }))
    group.length <- unlist(lapply(seq_len(max(group.name)), function(x) {
        return(length(which(group.name == x)))
    }))

    # RefSeqID
    RefSeqID <- trans.id[group.min]

    # targetName
    group.seq <- seqnames(A)
    targetName <- group.seq[group.min]
    targetName <- as.character(targetName)

    # strand
    group.strand <- strand(A)
    strand <- group.strand[group.min]
    strand <- as.character(strand)

    # targetStart
    targetStart <- lapply(group.index, function(x) {
        return(group.start[x])
    })

    # blockSizes
    blockSizes <- lapply(group.index, function(x) {
        return(group.width[x])
    })

    A.list <- vector2GRangesList(RefSeqID, targetName, strand, blockSizes, targetStart)

    return(A.list)
}

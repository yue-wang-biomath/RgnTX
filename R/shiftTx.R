#' Shift over transcripts
#' @export shiftTx
#'
#' @description Calculate positional shifting over transcript regions. This function accepts a feature set and outputs a region set from it. Each output region is from each input feature.
#'
#' @usage shiftTx(regions, start, width, direction, strand)
#'
#' @param regions A feature set following the format indicated in vignette section 3. Either to be \code{GRanges} or \code{GRangesList}.
#' @param start Starting positions. Each value represents a starting position in each input feature.
#' @param width Widths. Each value represents a width of each region to be picked from each feature.
#' @param direction Either to be character "left" or "right", which means the direction to which the starting position is shifting. The former means moving to the direction of 5' while the latter means moving to 3'.
#' @param strand The strand type of the transcripts. It receives "+" or "-".
#'
#' @return
#' A \code{Granges} object.
#'
#' @examples
#' # Take five transcripts.
#' # Extract the last 200 nt regions from their CDS part.
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' trans.id.pstv <- c("170", "782", "974", "1364", "1387")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' # download the CDS part of all transcriptome
#' cds.tx0 <- cdsBy(txdb, use.names = FALSE)
#'
#' # pick the CDS part of these five transcripts
#' cds.p <- cds.tx0[trans.id.pstv]
#'
#' width <- 200
#' start <- as.numeric(max(end(cds.p)))
#' R.cds.last200 <- shiftTx(cds.p, start = start, width = width, direction = 'left', strand = "+")

shiftTx = function(regions, start, width, direction, strand){
    if (is(regions, "GRanges")) {
        regions <- GRanges2GRangesList(regions)
    }

    if (!is(regions, "GRangesList")) {
        stop("regions should be either GRanges or GRangesList.")
    }
    disp <- data.frame(
        start = start,
        distance = width - 1,
        names = names(regions))

    regions.shift <- calculateShift(regions, disp, direction = direction, strand = strand)
    regions.return <- extractRegions(regions, regions.shift, strand = strand)
    return(regions.return)
}










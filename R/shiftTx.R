#' Shift over transcripts.
#' @export shiftTx
#'
#' @description Calculate positional shift over transcript regions.
#'
#' @usage shiftTx(regions, start, width, direction, strand)
#'
#' @param regions A feature set that follows the format indicated in vignette section 3. Either to be GRanges or GRangesList.
#' @param start A vector of starting positions, in which each value must be relative to each input feature.
#' @param width  A vector of integers, the width of each region to be picked from each feature.
#' @param direction The direction of displacement. It has options 'left' and 'right'. 'left' means shifting to 5' while 'right' to 3'.
#' @param strand The strand type of the transcripts. It receives '+' or '-'.
#'
#' @return
#' A Granges object.
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
#' start = as.numeric(max(end(cds.p)))
#' R.cds.last200 <- shiftTx(cds.p, start = start, width = width, direction = 'left', strand = "+")

shiftTx = function(regions, start, width, direction, strand){
    disp <- data.frame(
        start = start,
        distance = width - 1,
        names = names(regions))

    regions.shift <- calculateShift(regions, disp, direction = direction, strand = strand)
    regions.return <- extractRegions(regions, regions.shift, strand = strand)
    return(regions.return)
}










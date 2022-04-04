#' getStopCodon
#' @export getStopCodon
#' @description Get stop codon regions for input transcripts. This is an example of customPick function.
#'
#' @usage getStopCodon(trans_ids, txdb, ...)
#'
#' @param trans_ids A character object of transcript ids.
#' @param txdb A TxDb object.
#' @param ... any additional parameters needed.
#'
#' @return A \code{numeric} object.
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' trans.ids <- c("170", "782", "974", "1364", "1387")
#' RS2 <- getStopCodon(trans.ids, txdb)
getStopCodon <- function(trans_ids, txdb, ...) {
    # This example customPick function will get the stop codon regions of each transcript.

    dist <- 100
    trans.id <- as.character(trans_ids)
    cds.tx0 <- cdsBy(txdb, use.names = FALSE)
    threeUTR.tx0 <- threeUTRsByTranscript(txdb, use.names = FALSE)

    # The stop codon region is defined to be the first 100 nt regions
    # of 3'UTR and the last 100 nt regions of CDS in mRNAs.  The last
    # 100 nt regions of CDS.
    cds.names <- as.character(intersect(names(cds.tx0), trans.id))
    cds <- cds.tx0[cds.names]

    # separate p/n strand
    cds.strand <- data.frame(unique(strand(cds)))[, "value"]
    names.p <- cds.names[which(cds.strand == "+")]
    names.n <- cds.names[which(cds.strand == "-")]

    # p.strand
    if (length(names.p) != 0) {
        cds.p <- cds.tx0[names.p]
        cds.p.start <- as.numeric(max(end(cds.p)))
        disp.p.l <- data.frame(
            start = cds.p.start, distance = dist - 1,
            names = names.p
        )
        R.p.l <- calculateShift(cds.p, disp.p.l, direction = "left", strand = "+")
        R.p.l.exons <- extractRegions(regions_A = cds.p, R = R.p.l, strand = "+")
    } else {
        R.p.l.exons <- c()
    }

    # n.strand
    if (length(names.n) != 0) {
        cds.n <- cds.tx0[names.n]
        cds.n.start <- as.numeric(min(start(cds.n)))
        disp.n.l <- data.frame(
            start = cds.n.start, distance = dist - 1,
            names = names.n
        )
        R.n.l <- calculateShift(cds.n, disp.n.l, direction = "left", strand = "-")
        R.n.l.exons <- extractRegions(regions_A = cds.n, R = R.n.l, strand = "-")
    } else {
        R.n.l.exons <- c()
    }

    # The first 100 nt regions of 3'UTR
    threeUTR.names <- as.character(intersect(names(threeUTR.tx0), trans.id))
    threeUTR <- threeUTR.tx0[threeUTR.names]

    # separate p/n strand
    threeUTR.strand <- data.frame(unique(strand(threeUTR)))[, "value"]
    names.p <- threeUTR.names[which(threeUTR.strand == "+")]
    names.n <- threeUTR.names[which(threeUTR.strand == "-")]

    # p.strand
    if (length(names.p) != 0) {
        threeUTR.p <- threeUTR.tx0[names.p]
        threeUTR.p.start <- as.numeric(min(start(threeUTR.p)))
        disp.p.r <- data.frame(start = threeUTR.p.start, distance = dist -
            1, names = names.p)
        R.p.r <- calculateShift(threeUTR.p, disp.p.r,
            direction = "right",
            strand = "+"
        )
        R.p.r.exons <- extractRegions(
            regions_A = threeUTR.p, R = R.p.r,
            strand = "+"
        )
    } else {
        R.p.r.exons <- c()
    }

    # n.strand
    if (length(names.n) != 0) {
        threeUTR.n <- threeUTR.tx0[names.n]
        threeUTR.n.start <- as.numeric(max(end(threeUTR.n)))
        disp.n.r <- data.frame(start = threeUTR.n.start, distance = dist -
            1, names = names.n)
        R.n.r <- calculateShift(threeUTR.n, disp.n.r,
            direction = "right",
            strand = "-"
        )
        R.n.r.exons <- extractRegions(
            regions_A = threeUTR.n, R = R.n.r,
            strand = "-"
        )
    } else {
        R.n.r.exons <- c()
    }

    # rearrange
    R.stop.codon.p <- c(R.p.l.exons, R.p.r.exons)
    R.stop.codon.n <- c(R.n.l.exons, R.n.r.exons)

    # stop codon
    R.stop.codon.n$group <- R.stop.codon.n$group + max(R.stop.codon.p$group)
    if (length(R.stop.codon.p$group) != 0 && length(R.stop.codon.n$group) !=
        0) {
        R.stop.codon <- c(R.stop.codon.p, R.stop.codon.n)
    } else {
        if (length(R.stop.codon.p$group) != 0) {
            R.stop.codon <- R.stop.codon.p
        }
        if (length(R.stop.codon.n$group) != 0) {
            R.stop.codon <- R.stop.codon.n
        }
    }

    R.stop.codon <- GRanges2GRangesList(R.stop.codon)
    return(R.stop.codon)
}

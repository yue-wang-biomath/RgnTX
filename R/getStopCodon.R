#' getStopCodon
#' @export getStopCodon
#' @description Get stop codon regions for input transcripts. This is an example of customPick function.
#' @importFrom IRanges reduce
#' @importFrom GenomicRanges GRangesList
#' @usage getStopCodon(trans_ids, txdb, ...)
#'
#' @param trans_ids A character object containing transcript ids.
#' @param txdb A TxDb object.
#' @param ... Any additional parameters needed.
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
    trans.id <- unique(trans.id) ###
    cds.tx0 <- cdsBy(txdb, use.names = FALSE)
    threeUTR.tx0 <- threeUTRsByTranscript(txdb, use.names = FALSE)

    # The stop codon region is defined to be the first 100 nt regions
    # of 3'UTR and the last 100 nt regions of CDS in mRNAs.  The last
    # 100 nt regions of CDS.

    cds.names <- as.character(intersect(names(cds.tx0), trans.id))
    cds <- cds.tx0[cds.names]

    # separate p/n strand
    cds.strand <- data.frame(unique(strand(cds)))[, "value"]
    names.p <- cds.names[cds.strand == "+"]
    names.n <- cds.names[cds.strand == "-"]

    # p.strand
    if (length(names.p) != 0) {
        cds.p <- cds.tx0[names.p]
        cds.p.start <- as.numeric(max(end(cds.p)))
        R.p.l.exons <- shiftTx(cds.p, cds.p.start, dist, direction = "left", strand = "+")
        R.p.l.exons$group <- unlist(lapply(R.p.l.exons$transcriptsHits, function(x){return(which(trans.id == x))}))
    } else {
        R.p.l.exons <- c()
    }

    # n.strand
    if (length(names.n) != 0) {
        cds.n <- cds.tx0[names.n]
        cds.n.start <- as.numeric(min(start(cds.n)))
        R.n.l.exons <- shiftTx(cds.n, cds.n.start, dist, direction = "left", strand = "-")
        R.n.l.exons$group <- unlist(lapply(R.n.l.exons$transcriptsHits, function(x){return(which(trans.id == x))}))
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
        R.p.r.exons <- shiftTx(threeUTR.p, threeUTR.p.start, dist, direction = "right", strand = "+")
        R.p.r.exons$group <- unlist(lapply(R.p.r.exons$transcriptsHits, function(x){return(which(trans.id == x))}))

    } else {
        R.p.r.exons <- c()
    }

    # n.strand
    if (length(names.n) != 0) {
        threeUTR.n <- threeUTR.tx0[names.n]
        threeUTR.n.start <- as.numeric(max(end(threeUTR.n)))
        R.n.r.exons <- shiftTx(threeUTR.n, threeUTR.n.start, dist, direction = "right", strand = "-")
        R.n.r.exons$group <- unlist(lapply(R.n.r.exons$transcriptsHits, function(x){return(which(trans.id == x))}))

    } else {
        R.n.r.exons <- c()
    }

    # rearrange
    R.stop.codon <- c(R.p.l.exons, R.p.r.exons, R.n.l.exons, R.n.r.exons)
    groupall <- R.stop.codon$group
    groupindex <-unlist(lapply(unique(groupall), function(x){return(which(groupall == x))}))
    # stop codon
    R.stop.codon <- R.stop.codon[groupindex]
    groupname <- R.stop.codon$group
    R.stop.codon$group <- unlist(lapply(1:length(unique(groupname)),
                                        function(x){return(replicate(length(which(groupname == unique(groupname)[x])), x))}))

    R.stop.codon <- GRanges2GRangesList(R.stop.codon)
    R.stop.codon <- lapply(R.stop.codon, reduce)
    R.stop.codon <- GRangesList(R.stop.codon)
    return(R.stop.codon)
}

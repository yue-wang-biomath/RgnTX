#' Extract regions
#' @export extractRegions
#'
#' @description This function receives three arguments: the scope region set, the target region set and the type of strand. It returns a subset of target region set, which is the intersection of the target region set and the scope region set.
#'
#' @usage extractRegions(regions_A, R, strand = "+")
#'
#' @param regions_A The scope region set. A \code{GRangesList} object. The name of each list element should be the transcript id that it pertains to.
#' @param R The target region set. A \code{GRanges} object.
#' @param strand The strand type of the transcripts. It has options "+" and "-".
#'
#' @return A \code{GRangesList} object.
#'
#' @seealso \code{\link{calculateShift}}
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
#' disp.p.l <- data.frame(
#'     start = as.numeric(max(end(cds.p))),
#'     distance = width - 1,
#'     names = trans.id.pstv
#' )
#' R.p.l <- calculateShift(regions = cds.p, disp = disp.p.l, direction = "left", strand = "+")
#'
#' R.cds.last200 <- extractRegions(regions_A = cds.p, R = R.p.l, strand = "+")
extractRegions <- function(regions_A, R, strand = "+") {
    if (strand == "+") {
        start.pos <- start(R)
        end.pos <- end(R)
        regions.A <- regions_A
        A.end <- end(regions.A)
        A.start <- start(regions.A)

        R.start <- A.start[start.pos <= A.end]
        R.start <- R.start[R.start <= end.pos]
        R.start.frame <- data.frame(R.start)
        transcriptsHits <- as.character(R.start.frame[, "group_name"])

        R.end <- A.end[end.pos >= A.start]
        R.end <- R.end[R.end >= start.pos]
        R.end.frame <- data.frame(R.end)

        index.group <- unique(R.start.frame[, 1])
        R.start <- lapply(index.group, function(x) {
            regs <- R.start.frame[R.start.frame[, 1] == x, ]
            regs[1, 3] <- start.pos[x]
            return(regs[, 3])
        })
        R.end <- lapply(index.group, function(x) {
            regs <- R.end.frame[R.end.frame[, 1] == x, ]
            regs[nrow(regs), 3] <- end.pos[x]
            return(regs[, 3])
        })
        R.start.frame[, 3] <- unlist(R.start)
        R.end.frame[, 3] <- unlist(R.end)

        regions.R <- regions.A[as.character(R.start.frame[, 2])]
        R <- GRanges(
            seqnames = Rle(unlist(unique(seqnames(regions.R)))),
            IRanges(start = R.start.frame[, 3], end = R.end.frame[, 3]),
            strand = "+", transcriptsHits = transcriptsHits, group = R.start.frame[
                ,
                "group"
            ]
        )
    }
    if (strand == "-") {
        regions.A <- regions_A
        start.pos <- start(R)
        end.pos <- end(R)
        A.end <- end(regions.A)
        A.start <- start(regions.A)

        R.start <- A.start[start.pos <= A.end]
        R.start <- R.start[R.start <= end.pos]
        R.start.frame <- data.frame(R.start)
        transcriptsHits <- as.character(R.start.frame[, "group_name"])

        R.end <- A.end[end.pos >= A.start]
        R.end <- R.end[R.end >= start.pos]
        R.end.frame <- data.frame(R.end)

        index.group <- unique(R.start.frame[, 1])
        R.start <- lapply(index.group, function(x) {
            regs <- R.start.frame[R.start.frame[, 1] == x, ]
            regs[nrow(regs), 3] <- start.pos[x]
            return(regs[, 3])
        })
        R.end <- lapply(index.group, function(x) {
            regs <- R.end.frame[R.end.frame[, 1] == x, ]
            regs[1, 3] <- end.pos[x]
            return(regs[, 3])
        })
        R.start.frame[, 3] <- unlist(R.start)
        R.end.frame[, 3] <- unlist(R.end)

        regions.R <- regions.A[as.character(R.start.frame[, 2])]
        R <- GRanges(
            seqnames = Rle(unlist(unique(seqnames(regions.R)))),
            IRanges(start = R.start.frame[, 3], end = R.end.frame[, 3]),
            strand = "-", transcriptsHits = transcriptsHits, group = R.start.frame[
                ,
                "group"
            ]
        )
    }
    return(R)
}

#' Calculate positional shifting over transcriptome
#' @description This function can calculate shifting over 'disconnected' regions. For example, due to alternative splicing, an mRNA is represented by several disconnected regions according to the annotation information. And this function can calculate position displacement over such disconnected regions. It can automatically skip the 'blank' between the regions when calculating position shifting.
#' @usage calculateShift(regions, disp, direction = 'right', strand = '+')
#'
#' @param regions The region space where the displacement is calculated, which should be a GRangesList object.
#' @param disp The region space where the displacement is calculated, which should be a GRangesList object.
#' @param direction The direction of displacement. It has options 'left' and 'right'.
#' @param strand The strand type of the transcripts. It receives '+' or '-'.
#'
#' @export calculateShift
#' @importFrom IRanges IntegerList IRanges LogicalList
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges end makeGRangesListFromFeatureFragments
#' @importFrom S4Vectors Rle
#' @import GenomicAlignments
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' @return A \code{GRanges} object.
#'
#' @seealso \code{\link{extractRegions}}
#'
#' @examples
#' # Take five transcripts.
#' # Extract the last 200 nt regions from their CDS part.
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' trans.id.pstv <- c("170", "782", "974", "1364", "1387")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' # Download the CDS part of all transcriptome
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
#' R.p.l <- calculateShift(
#'     regions = cds.p, disp = disp.p.l,
#'     direction = "left", strand = "+"
#' )
calculateShift <-
    function(regions, disp, direction = "right", strand = "+") {
        if (direction == "right" && strand == "+") {
            start.pos <- disp[, 1]
            displacement <- disp[, 2]
            regions.A <- regions
            A.end <- end(regions.A)
            A.start <- start(regions.A)

            index.reg <- lapply(start.pos <= A.end, function(x) {
                return(which(x))
            })

            reg1.width <- start.pos - min(A.start[index.reg])
            reg2.width <- reg1.width + displacement
            width.cumsum <- cumsum(width(regions.A)[index.reg])
            index.reg.1 <- (reg2.width >= width.cumsum)

            index.reg.1 <- LogicalList(lapply(index.reg.1, function(x) {
                x <- c(TRUE, x)
                x <- x[seq_len(length(x) - 1)]
                return(x)
            }))

            reg3.width <- max(width.cumsum[index.reg.1])
            reg4.width <- reg3.width - reg2.width
            end.pos <- max(A.end[index.reg][index.reg.1]) - reg4.width + 1
            end.pos[reg2.width > max(width.cumsum)] <- max(A.end[index.reg][index.reg.1])[reg2.width > max(width.cumsum)]

            R <- GRanges(
                seqnames = Rle(unlist(unique(seqnames(regions.A)))),
                IRanges(
                    start = start.pos,
                    end = end.pos
                ),
                strand = "+",
                transcriptsHits = as.character(disp[, "names"]),
                group = seq_len(nrow(disp))
            )
        }
        if (direction == "right" && strand == "-") {
            # run
            start.pos <- disp[, 1]
            displacement <- disp[, 2]
            regions.A <- regions
            A.end <- end(regions.A)
            A.start <- start(regions.A)

            index.reg <- lapply(start.pos >= A.start, function(x) {
                return(which(x))
            })
            reg1.width <- max(A.end[index.reg]) - start.pos
            reg2.width <- reg1.width + displacement
            width.cumsum <- cumsum(width(regions.A)[index.reg])

            index.reg.1 <- (reg2.width >= width.cumsum)
            index.reg.1 <- LogicalList(lapply(index.reg.1, function(x) {
                x <- c(TRUE, x)
                x <- x[seq_len(length(x) - 1)]
                return(x)
            }))

            reg3.width <- max(width.cumsum[index.reg.1])
            reg4.width <- reg3.width - reg2.width
            end.pos <- min(A.start[index.reg][index.reg.1]) + reg4.width - 1
            end.pos[reg2.width > max(width.cumsum)] <- min(A.start[index.reg][index.reg.1])[reg2.width > max(width.cumsum)]

            R <- GRanges(
                seqnames = Rle(unlist(unique(seqnames(regions.A)))),
                IRanges(
                    start = end.pos,
                    end = start.pos
                ),
                strand = "-",
                transcriptsHits = as.character(disp[, "names"]),
                group = seq_len(nrow(disp))
            )
        }
        if (direction == "left" && strand == "+") {
            start.pos <- disp[, 1]
            displacement <- disp[, 2]
            regions.A <- regions
            A.end <- end(regions.A)
            A.start <- start(regions.A)

            index.reg <- lapply(start.pos >= A.start, function(x) {
                return(which(x))
            })
            reg1.width <- max(A.end[index.reg]) - start.pos
            reg2.width <- displacement + reg1.width
            width.cumsum <- cumsum(width(regions.A)[lapply(index.reg, rev)])
            width.cumsum <- IntegerList(lapply(width.cumsum, rev))
            index.reg.1 <- reg2.width >= width.cumsum

            index.reg.1 <- LogicalList(lapply(index.reg.1, function(x) {
                x <- rev(x)
                x <- c(TRUE, x)
                x <- x[seq_len(length(x) - 1)]
                return(rev(x))
            }))

            reg3.width <- max(width.cumsum[index.reg.1])
            reg4.width <- reg3.width - reg2.width
            end.pos <- min(A.start[index.reg][index.reg.1]) + reg4.width - 1
            end.pos[reg2.width > max(width.cumsum)] <- min(A.start[index.reg][index.reg.1])[reg2.width > max(width.cumsum)]


            R <- GRanges(
                seqnames = Rle(unlist(unique(seqnames(regions.A)))),
                IRanges(
                    start = end.pos,
                    end = start.pos
                ),
                strand = "+",
                transcriptsHits = as.character(disp[, "names"]),
                group = seq_len(nrow(disp))
            )
        }
        if (direction == "left" && strand == "-") {
            # run
            start.pos <- disp[, 1]
            displacement <- disp[, 2]
            regions.A <- regions
            A.end <- end(regions.A)
            A.start <- start(regions.A)

            index.reg <- lapply(start.pos <= A.end, function(x) {
                return(which(x))
            })
            reg1.width <- start.pos - min(A.start[index.reg])
            reg2.width <- reg1.width + displacement
            width.cumsum <- cumsum(width(regions.A)[lapply(index.reg, function(x) {
                return(rev(x))
            })])
            width.cumsum <- IntegerList(lapply(width.cumsum, function(x) {
                return(rev(x))
            }))
            index.reg.1 <- reg2.width >= width.cumsum
            index.reg.1 <- LogicalList(lapply(index.reg.1, function(x) {
                x <- rev(x)
                x <- c(TRUE, x)
                x <- x[seq_len(length(x) - 1)]
                return(rev(x))
            }))
            reg3.width <- max(width.cumsum[index.reg.1])
            reg4.width <- reg3.width - reg2.width
            end.pos <- max(A.end[index.reg][index.reg.1]) - reg4.width + 1

            end.pos[reg2.width > max(width.cumsum)] <- max(A.end)[reg2.width > max(width.cumsum)]

            R <- GRanges(
                seqnames = Rle(unlist(unique(seqnames(regions.A)))),
                IRanges(
                    start = start.pos,
                    end = end.pos
                ),
                strand = "-",
                transcriptsHits = as.character(disp[, "names"]),
                group = seq_len(nrow(disp))
            )
        }

        return(R)
    }

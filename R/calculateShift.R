#' Calculate positional shifting over transcriptome
#' @description The first step of calculating positional shift over transcriptome regions.
#' @usage calculateShift(regions, disp, direction = "right", strand = "+")
#'
#' @param regions A feature set, which should be a GRangesList object.
#' @param disp A data frame object. It should have three columns, which are \code{start}: starting positions. Each value represents a starting position in each input feature;
#' \code{width}: widths. Each value represents a width of each region to be picked from each feature; \code{names}: corresponding transcript ids.
#' @param direction Either to be character "left" or "right", which means the direction to which the starting position is shifting. The former means moving to the direction of 5' while the latter means moving to 3'.
#' @param strand Either to be "+" or "-".
#'
#' @export calculateShift
#' @importFrom IRanges IntegerList IRanges LogicalList shift
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges end makeGRangesListFromFeatureFragments
#' @importFrom S4Vectors Rle
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
            rightmost <- max(A.end[index.reg][index.reg.1])
            end.pos <- rightmost - reg4.width + 1
            end.pos[reg2.width > max(width.cumsum)] <-rightmost[reg2.width > max(width.cumsum)]
            end.pos[rightmost < end.pos] <- rightmost[rightmost < end.pos]

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
            rightmost <- min(A.start[index.reg][index.reg.1])
            end.pos <-  rightmost + reg4.width - 1
            end.pos[reg2.width > max(width.cumsum)] <- rightmost[reg2.width > max(width.cumsum)]
            end.pos[rightmost > end.pos] <- rightmost[rightmost > end.pos]

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
            leftmost <- min(A.start[index.reg][index.reg.1])
            end.pos <- leftmost + reg4.width - 1
            end.pos[reg2.width > max(width.cumsum)] <- leftmost[reg2.width > max(width.cumsum)]
            end.pos[leftmost > end.pos] <- leftmost[leftmost > end.pos]


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
            leftmost <- max(A.end[index.reg][index.reg.1])
            end.pos <- leftmost - reg4.width + 1
            end.pos[reg2.width > max(width.cumsum)] <- leftmost[reg2.width > max(width.cumsum)]
            end.pos[leftmost < end.pos] <- leftmost[leftmost < end.pos]


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

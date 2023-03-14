#' Shift regions over transcripts
#' @export shiftExonTx
#' @importFrom GenomicRanges GRangesList
#' @description Shift regions over transcripts (exons). This function accepts a feature set and outputs a region set from it. Each output region is from each input feature.
#'
#' @usage shiftExonTx(regions, start, width)
#'
#' @param regions A region set either to be \code{GRanges} or \code{GRangesList}, following the format indicated in vignette section 3.
#' @param start Starting positions. The function will return a region set shifting from these starting positions.
#' @param width Widths. The widths of returned regions. It should be positive/negative if users would like to move the starting positions to the right (3')/left (5').
#'
#' @return
#' A \code{GRangesList} object.
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' exons.tx0 <- exonsBy(txdb)
#' #Shift regions over the following transcripts
#' trans.ids <- c("170", "1387", "4113", "10715")
#' regions <- exons.tx0[trans.ids]
#' start <- c(3624255, 55158197, 881641, 15694195)
#' width <- c(200, -200, 200, -200)
#' shifted_regions <- shiftExonTx(regions, start, width)
shiftExonTx = function(regions, start, width){
    if (is(regions, "GRanges")) {
        regions <- GRanges2GRangesList(regions)
    }

    if (!is(regions, "GRangesList")) {
        stop("regions should be either GRanges or GRangesList.")
    }

    direction <- replicate(length(regions), 'right')
    direction[width < 0] <- 'left'
    dataframe.r <- data.frame(regions)
    strand <- dataframe.r[, 'strand']
    strand <- strand[unlist(lapply(unique(dataframe.r[,1]), function(x){return(min(which(dataframe.r[, 1] == x)))}))]

    index_posv_right <- intersect(which(direction == 'right'), which(strand == '+'))
    index_posv_left <- intersect(which(direction == 'left'), which(strand == '+'))
    index_ngtv_right <- intersect(which(direction == 'right'), which(strand == '-'))
    index_ngtv_left <- intersect(which(direction == 'left'), which(strand == '-'))

    if(length(index_posv_right)!= 0){
        regions.return_posv_right <- shiftTx(regions[index_posv_right], start[index_posv_right], width[index_posv_right], 'right', '+')
        group1 <- regions.return_posv_right$group
        group1 <- unlist(lapply(1:length(unique(group1)),
                         function(x){return(replicate(length(which(group1 ==unique(group1)[x])), x))}))
        regions.return_posv_right$group <- group1
        regions.return_posv_right <- GRanges2GRangesList(regions.return_posv_right)
    }else{
        regions.return_posv_right <- GRangesList()
    }


    if(length(index_posv_left)!= 0){
        regions.return_posv_left <- shiftTx(regions[index_posv_left], start[index_posv_left], -width[index_posv_left], 'left', '+')
        group2 <- regions.return_posv_left$group
        group2 <- unlist(lapply(1:length(unique(group2)),
                                function(x){return(replicate(length(which(group2 ==unique(group2)[x])), x))}))
        regions.return_posv_left$group <- group2
        regions.return_posv_left <- GRanges2GRangesList(regions.return_posv_left)
    }else{
        regions.return_posv_left <- GRangesList()
    }


    if(length(index_ngtv_right)!= 0){
        regions.return_ngtv_right <- shiftTx(regions[index_ngtv_right], start[index_ngtv_right], width[index_ngtv_right], 'right', '-')
        group3 <- regions.return_ngtv_right$group
        group3 <- unlist(lapply(1:length(unique(group3)),
                                function(x){return(replicate(length(which(group3 ==unique(group3)[x])), x))}))
        regions.return_ngtv_right$group <- group3
        regions.return_ngtv_right <- GRanges2GRangesList(regions.return_ngtv_right)
    }else{
        regions.return_ngtv_right <- GRangesList()
    }

    if(length(index_ngtv_left)!= 0){
        regions.return_ngtv_left <- shiftTx(regions[index_ngtv_left], start[index_ngtv_left], -width[index_ngtv_left], 'left', '-')
        group4 <- regions.return_ngtv_left$group
        group4 <- unlist(lapply(1:length(unique(group4)),
                                function(x){return(replicate(length(which(group4 ==unique(group4)[x])), x))}))
        regions.return_ngtv_left$group <- group4
        regions.return_ngtv_left <- GRanges2GRangesList(regions.return_ngtv_left)
    }else{
        regions.return_ngtv_left <- GRangesList()
    }

    regions.return <- c(regions.return_posv_right, regions.return_posv_left, regions.return_ngtv_right, regions.return_ngtv_left)

    return(regions.return)
}
